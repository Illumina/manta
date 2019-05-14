//
// Manta - Structural Variant and Indel Caller
// Copyright (c) 2013-2019 Illumina, Inc.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//

/// \file
/// \author Chris Saunders
///

#include "GenerateSVCandidates.hpp"

#include <iostream>
#include <string>

#include "EdgeRetrieverBin.hpp"
#include "EdgeRetrieverLocus.hpp"
#include "GSCOptions.hpp"
#include "SVCandidateProcessor.hpp"
#include "SVEvidenceWriter.hpp"
#include "SVFinder.hpp"
#include "blt_util/log.hpp"
#include "common/Exceptions.hpp"
#include "manta/BamStreamerUtils.hpp"
#include "manta/MultiJunctionUtil.hpp"
#include "manta/SVCandidateUtil.hpp"

// A bus error on PPC64 (DRAGEN-1789) seems to be occurring in the boost lockfree queue logic
// of the boost version of CTPL, although this was not traced down to a specific error mechanism.
// To avoid this error under the assumption that there is a problem with boost CTPL, switch to
// the STL version of CTPL instead:
#include "ctpl_stl.h"

//#define DEBUG_GSV

/// Provide additional edge details, intended for attachment to an in-flight exception:
static void dumpEdgeInfo(const EdgeInfo& edge, const SVLocusSet& set, std::ostream& os)
{
  const bam_header_info& bamHeader(set.getBamHeader());
  os << edge;
  const auto& node1(set.getLocus(edge.locusIndex).getNode(edge.nodeIndex1));
  os << "\tnode1:" << node1;
  os << "\tnode1:";
  summarizeGenomeInterval(bamHeader, node1.getInterval(), os);
  os << "\n";

  if (!edge.isSelfEdge()) {
    const auto& node2(set.getLocus(edge.locusIndex).getNode(edge.nodeIndex2));
    os << "\tnode2:" << node2;
    os << "\tnode2:";
    summarizeGenomeInterval(bamHeader, node2.getInterval(), os);
    os << "\n";
  }
}

/// we can either traverse all edges in a single locus (disjoint subgraph) of the graph
/// OR
/// traverse all edges in one "bin" -- that is, one out of binCount subsets of the total
/// graph edges. Each bin is designed to be of roughly equal size in terms of total
/// anticipated workload, so that we have good parallel processing performance.
///
static EdgeRetriever* edgeRFactory(const SVLocusSet& set, const EdgeOptions& opt)
{
  if (opt.isLocusIndex) {
    return (new EdgeRetrieverLocus(set, opt.graphNodeMaxEdgeCount, opt.locusOpt));
  } else {
    return (new EdgeRetrieverBin(set, opt.graphNodeMaxEdgeCount, opt.binCount, opt.binIndex));
  }
}

/// TODO temporarily shoved here, needs a better home:
static void multiJunctionFilterGroupCandidateSV(
    const GSCOptions&                      opt,
    const EdgeInfo&                        edge,
    const std::vector<SVCandidate>&        svs,
    GSCEdgeStatsManager&                   edgeStatMan,
    std::vector<SVMultiJunctionCandidate>& mjSVs)
{
  unsigned mjComplexCount(0);
  unsigned mjSpanningFilterCount(0);
  findMultiJunctionCandidates(
      svs, opt.minCandidateSpanningCount, opt.isRNA, mjComplexCount, mjSpanningFilterCount, mjSVs);
  edgeStatMan.updateMJFilter(edge, mjComplexCount, mjSpanningFilterCount);

  if (opt.isVerbose) {
    unsigned junctionCount(mjSVs.size());
    unsigned candidateCount(0);
    for (const SVMultiJunctionCandidate& mj : mjSVs) {
      candidateCount += mj.junction.size();
    }
    log_os << __FUNCTION__ << ": Low-resolution candidate filtration complete. "
           << "candidates: " << candidateCount << " "
           << "junctions: " << junctionCount << " "
           << "complex: " << mjComplexCount << " "
           << "spanningfilt: " << mjSpanningFilterCount << "\n";
  }
#ifdef DEBUG_GSV
  log_os << __FUNCTION__ << ": final candidate list";
  const unsigned junctionCount(mjSVs.size());
  for (unsigned junctionIndex(0); junctionIndex < junctionCount; ++junctionIndex) {
    const auto&    mj(mjSVs[junctionIndex]);
    const unsigned junctionCandCount(mj.junction.size());
    log_os << __FUNCTION__ << ": JUNCTION " << junctionIndex << " with " << junctionCandCount
           << " candidates\n";
    for (unsigned junctionCandIndex(0); junctionCandIndex < junctionCandCount; ++junctionCandIndex) {
      log_os << __FUNCTION__ << ":  JUNCTION " << junctionIndex << " Candidate " << junctionCandIndex << " "
             << mj.junction[junctionCandIndex] << "\n";
    }
  }
#endif
}

namespace {

/// Mutable data shared by all worker threads:
std::atomic<bool> isWorkerThreadException(false);

}  // namespace

/// These are thread-local data that we initialize before the thread pool to reduce total initialization
/// costs, specifically we want to avoid initializing them for every edge (in principal, this has not been
/// benchmarked).
struct EdgeThreadLocalData {
  std::shared_ptr<EdgeRuntimeTracker>   edgeTrackerPtr;
  GSCEdgeStatsManager                   edgeStatMan;
  std::unique_ptr<SVFinder>             svFindPtr;
  std::unique_ptr<SVCandidateProcessor> svProcessorPtr;
  std::unique_ptr<SVEvidenceWriter>     svEvidenceWriterPtr;
  SVCandidateSetData                    svData;
  std::vector<SVCandidate>              svs;
  std::vector<SVMultiJunctionCandidate> mjSVs;
};

/// Process a single edge on one thread:
static void processEdge(
    int                               threadId,
    const GSCOptions&                 opt,
    const SVLocusSet&                 cset,
    std::vector<EdgeThreadLocalData>& edgeDataPool,
    const EdgeInfo                    edge)
{
  if (isWorkerThreadException.load()) return;

  EdgeThreadLocalData& edgeData(edgeDataPool[threadId]);

  try {
    edgeData.edgeTrackerPtr->start();

    if (opt.isVerbose) {
      log_os << __FUNCTION__ << ": starting analysis of edge: ";
      dumpEdgeInfo(edge, cset, log_os);
    }

    // find number, type and breakend range (or better: breakend distro) of SVs on this edge:
    edgeData.svFindPtr->findCandidateSV(cset, edge, edgeData.svData, edgeData.svs);

    // filter long-range junctions outside of the candidate finder so that we can evaluate
    // junctions which are part of a larger event (like a reciprocal translocation)
    multiJunctionFilterGroupCandidateSV(opt, edge, edgeData.svs, edgeData.edgeStatMan, edgeData.mjSVs);

    // assemble, score and output SVs
    edgeData.svProcessorPtr->evaluateCandidates(edge, edgeData.mjSVs, edgeData.svData);
  } catch (illumina::common::ExceptionData& e) {
    isWorkerThreadException = true;
    std::ostringstream oss;
    oss << "Exception caught in thread " << threadId << " while processing graph edge: ";
    dumpEdgeInfo(edge, cset, oss);
    e << boost::error_info<struct current_edge_info, std::string>(oss.str());
    throw;
  } catch (std::exception& e) {
    isWorkerThreadException = true;

    // Convert std::exception to boost exception so that we can add on more context:
    illumina::common::GeneralException ge(e.what());

    std::ostringstream oss;
    oss << "Exception caught in thread " << threadId << " while processing graph edge: ";
    dumpEdgeInfo(edge, cset, oss);
    ge << boost::error_info<struct current_edge_info, std::string>(oss.str());
    BOOST_THROW_EXCEPTION(ge);
  } catch (...) {
    log_os << "Unknown exception caught while processing graph edge: ";
    dumpEdgeInfo(edge, cset, log_os);
    throw;
  }

  edgeData.edgeTrackerPtr->stop(edge);
  if (opt.isVerbose) {
    log_os << __FUNCTION__ << ": Time to process last edge: ";
    edgeData.edgeTrackerPtr->getLastEdgeTime().reportSec(log_os);
    log_os << "\n";
  }

  edgeData.edgeStatMan.updateScoredEdgeTime(edge, *(edgeData.edgeTrackerPtr));
}

static void runGSC(const GSCOptions& opt, const char* progName, const char* progVersion)
{
  const SVLocusScanner readScanner(
      opt.scanOpt, opt.statsFilename, opt.alignFileOpt.alignmentFilenames, !opt.isUnstrandedRNA);

  static const bool      isSkipLocusSetIndexCreation(true);
  const SVLocusSet       cset(opt.graphFilename.c_str(), isSkipLocusSetIndexCreation);
  const bam_header_info& bamHeader(cset.getBamHeader());

  if (opt.isVerbose) {
    log_os << __FUNCTION__ << ": " << bamHeader << "\n";
  }

  const SVWriter svWriter(opt, bamHeader, progName, progVersion);
  auto           svEvidenceWriterSharedData(std::make_shared<SVEvidenceWriterSharedData>(opt));

  // Initialize all thread-local edge data:
  std::shared_ptr<SynchronizedOutputStream> edgeTrackerStreamPtr;
  if (!opt.edgeRuntimeFilename.empty()) {
    edgeTrackerStreamPtr.reset(new SynchronizedOutputStream(opt.edgeRuntimeFilename));
  }

  std::vector<EdgeThreadLocalData> edgeDataPool(opt.workerThreadCount);
  for (auto& edgeData : edgeDataPool) {
    edgeData.edgeTrackerPtr.reset(new EdgeRuntimeTracker(edgeTrackerStreamPtr));
    edgeData.svFindPtr.reset(new SVFinder(
        opt,
        readScanner,
        bamHeader,
        cset.getAllSampleReadCounts(),
        edgeData.edgeTrackerPtr,
        edgeData.edgeStatMan));
    edgeData.svProcessorPtr.reset(new SVCandidateProcessor(
        opt,
        readScanner,
        cset,
        svWriter,
        svEvidenceWriterSharedData,
        edgeData.edgeTrackerPtr,
        edgeData.edgeStatMan));
  }

  ctpl::thread_pool pool(opt.workerThreadCount);

  // Iterate through graph edges:
  std::unique_ptr<EdgeRetriever> edgerPtr(edgeRFactory(cset, opt.edgeOpt));
  EdgeRetriever&                 edger(*edgerPtr);

  // Although processEdge doesn't return anything now, the future<void> provides a simple way for worker
  // thread
  // exceptions to propogate down to this thread:
  std::vector<std::future<void>> edgeReturnValues;

  while (edger.next()) {
    edgeReturnValues.push_back(
        pool.push(processEdge, std::cref(opt), std::cref(cset), std::ref(edgeDataPool), edger.getEdge()));
  }

  pool.stop(true);

  // This is sufficient to rethrow any worker thread exceptions.
  //
  // Note that we don't catch them here, just throw down to the bottom- level handler for
  // GenerateSVCandidates:
  for (auto& edgeReturnValue : edgeReturnValues) {
    edgeReturnValue.get();
  }

  GSCEdgeStats mergedStats;
  for (auto& edgeData : edgeDataPool) {
    mergedStats.merge(edgeData.edgeStatMan.returnStats());
  }

  if (!opt.edgeStatsFilename.empty()) {
    mergedStats.save(opt.edgeStatsFilename.c_str());
  }
  if (!opt.edgeStatsReportFilename.empty()) {
    mergedStats.report(opt.edgeStatsReportFilename.c_str());
  }
}

void GenerateSVCandidates::runInternal(int argc, char* argv[]) const
{
  GSCOptions opt;

  parseGSCOptions(*this, argc, argv, opt);
#ifdef DEBUG_GSV
  opt.isVerbose = true;
#endif
  runGSC(opt, name(), version());
}
