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

#include "GenerateSVCandidates.hh"
#include "EdgeRetrieverBin.hh"
#include "EdgeRetrieverLocus.hh"
#include "GSCOptions.hh"
#include "SVCandidateProcessor.hh"
#include "SVFinder.hh"
#include "SVSupports.hh"

#include "blt_util/log.hh"
#include "common/Exceptions.hh"
#include "manta/BamStreamerUtils.hh"
#include "manta/MultiJunctionUtil.hh"
#include "manta/SVCandidateUtil.hh"

#include "ctpl.h"

#include <iostream>
#include <string>

//#define DEBUG_GSV


/// provide additional edge details, intended for attachment to an in-flight exception:
static
void
dumpEdgeInfo(
    const EdgeInfo& edge,
    const SVLocusSet& set,
    std::ostream& os)
{
    const bam_header_info& bamHeader(set.getBamHeader());
    os << edge;
    const auto& node1(set.getLocus(edge.locusIndex).getNode(edge.nodeIndex1));
    os << "\tnode1:" << node1;
    os << "\tnode1:";
    summarizeGenomeInterval(bamHeader, node1.getInterval(), os);
    os << "\n";

    if (! edge.isSelfEdge())
    {
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
static
EdgeRetriever*
edgeRFactory(
    const SVLocusSet& set,
    const EdgeOptions& opt)
{
    if (opt.isLocusIndex)
    {
        return (new EdgeRetrieverLocus(set, opt.graphNodeMaxEdgeCount, opt.locusOpt));
    }
    else
    {
        return (new EdgeRetrieverBin(set, opt.graphNodeMaxEdgeCount, opt.binCount, opt.binIndex));
    }
}



/// TODO temporarily shoved here, needs a better home:
static
void
multiJunctionFilterGroupCandidateSV(
    const GSCOptions& opt,
    const EdgeInfo& edge,
    const std::vector<SVCandidate>& svs,
    GSCEdgeStatsManager& edgeStatMan,
    std::vector<SVMultiJunctionCandidate>& mjSVs)
{
    unsigned mjComplexCount(0);
    unsigned mjSpanningFilterCount(0);
    findMultiJunctionCandidates(svs, opt.minCandidateSpanningCount, opt.isRNA, mjComplexCount, mjSpanningFilterCount, mjSVs);
    edgeStatMan.updateMJFilter(edge, mjComplexCount, mjSpanningFilterCount);

    if (opt.isVerbose)
    {
        unsigned junctionCount(mjSVs.size());
        unsigned candidateCount(0);
        for (const SVMultiJunctionCandidate& mj : mjSVs)
        {
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
    for (unsigned junctionIndex(0); junctionIndex< junctionCount; ++junctionIndex)
    {
        const auto& mj(mjSVs[junctionIndex]);
        const unsigned junctionCandCount(mj.junction.size());
        log_os << __FUNCTION__ << ": JUNCTION " << junctionIndex << " with " << junctionCandCount << " candidates\n";
        for (unsigned junctionCandIndex(0); junctionCandIndex< junctionCandCount; ++junctionCandIndex)
        {
            log_os << __FUNCTION__ << ":  JUNCTION " << junctionIndex << " Candidate "
                   << junctionCandIndex << " " << mj.junction[junctionCandIndex] << "\n";
        }
    }
#endif
}



/// These are thread-local data that we initialize before the thread pool to reduce total initialization costs,
/// specifically we want to avoid initializing them for every edge (in principal, this has not been benchmarked).
struct EdgeThreadLocalData {
    std::shared_ptr<EdgeRuntimeTracker> edgeTrackerPtr;
    GSCEdgeStatsManager edgeStatMan;
    std::unique_ptr<SVFinder> svFindPtr;
    SVCandidateSetData svData;
    std::vector<SVCandidate> svs;
    std::vector<SVMultiJunctionCandidate> mjSVs;
};



/// Process a single edge on one thread:
GSCEdgeStats
processEdge(
    int threadId,
    const GSCOptions& opt,
    const SVLocusSet& cset,
    std::vector<EdgeThreadLocalData>& edgeDataPool,
    const EdgeInfo edge)
{
    std::cerr << "Hello from thread " << threadId;
    EdgeThreadLocalData& edgeData(edgeDataPool[threadId]);

    try
    {
        edgeData.edgeTrackerPtr->start();

        if (opt.isVerbose)
        {
            log_os << __FUNCTION__ << ": starting analysis of edge: ";
            dumpEdgeInfo(edge,cset,log_os);
        }

        // find number, type and breakend range (or better: breakend distro) of SVs on this edge:
        edgeData.svFindPtr->findCandidateSV(cset, edge, edgeData.svData, edgeData.svs);

        // filter long-range junctions outside of the candidate finder so that we can evaluate
        // junctions which are part of a larger event (like a reciprocal translocation)
        multiJunctionFilterGroupCandidateSV(opt, edge, edgeData.svs,  edgeData.edgeStatMan, edgeData.mjSVs);

#ifdef RABBIT
        SupportSamples svSupports;
        svSupports.supportSamples.resize(sampleSize);

        // assemble, score and output SVs
        svProcessor.evaluateCandidates(edge, mjSVs, svData, svSupports);

        // write supporting reads into bam files
        if (isGenerateSupportBam)
        {
            for (unsigned idx(0); idx<sampleSize; ++idx)
            {
                writeSupportBam(origBamStreamPtrs[idx],
                                svSupports.supportSamples[idx],
                                supportBamDumperPtrs[idx]);
            }
        }
#endif
    }
    catch (illumina::common::ExceptionData& e)
    {
        std::ostringstream oss;
        oss << "Exception caught while processing graph edge: ";
        dumpEdgeInfo(edge,cset,oss);
        e << boost::error_info<struct current_edge_info,std::string>(oss.str());
        throw;
    }
    catch (...)
    {
        log_os << "Exception caught while processing graph edge: ";
        dumpEdgeInfo(edge,cset,log_os);
        throw;
    }

    edgeData.edgeTrackerPtr->stop(edge);
    if (opt.isVerbose)
    {
        log_os << __FUNCTION__ << ": Time to process last edge: ";
        edgeData.edgeTrackerPtr->getLastEdgeTime().reportSec(log_os);
        log_os << "\n";
    }

    edgeData.edgeStatMan.updateScoredEdgeTime(edge, *(edgeData.edgeTrackerPtr));
    return edgeData.edgeStatMan.returnStats();
}

static
void
runGSC(
    const GSCOptions& opt,
    const char* progName,
    const char* progVersion)
{
    auto edgeTrackerStreamPtr(std::make_shared<SynchronizedOutputStream>(opt.edgeRuntimeFilename));

    const SVLocusScanner readScanner(opt.scanOpt, opt.statsFilename, opt.alignFileOpt.alignmentFilenames, !opt.isUnstrandedRNA);

    static const bool isSkipLocusSetIndexCreation(true);
    const SVLocusSet cset(opt.graphFilename.c_str(), isSkipLocusSetIndexCreation);
    const bam_header_info& bamHeader(cset.getBamHeader());
#ifdef RABBIT
    SVCandidateProcessor svProcessor(opt, readScanner, progName, progVersion, cset, edgeTracker, edgeStatMan);
#endif

    std::unique_ptr<EdgeRetriever> edgerPtr(edgeRFactory(cset, opt.edgeOpt));
    EdgeRetriever& edger(*edgerPtr);



    const unsigned sampleSize(opt.alignFileOpt.alignmentFilenames.size());
    std::vector<bam_streamer_ptr> origBamStreamPtrs;
    std::vector<bam_dumper_ptr> supportBamDumperPtrs;

    const bool isGenerateSupportBam(opt.supportBamStub.size() > 0);
    if (isGenerateSupportBam)
    {
        openBamStreams(opt.referenceFilename, opt.alignFileOpt.alignmentFilenames, origBamStreamPtrs);

        assert(origBamStreamPtrs.size() == sampleSize);
        for (unsigned sampleIndex(0); sampleIndex<sampleSize; ++sampleIndex)
        {
            std::string supportBamName(opt.supportBamStub
                                       + ".bam_" + std::to_string(sampleIndex)
                                       + ".bam");
            const bam_hdr_t& header(origBamStreamPtrs[sampleIndex]->get_header());
            bam_dumper_ptr bamDumperPtr(new bam_dumper(supportBamName.c_str(), header));
            supportBamDumperPtrs.push_back(bamDumperPtr);
        }
    }

    if (opt.isVerbose)
    {
        log_os << __FUNCTION__ << ": " << bamHeader << "\n";
    }

    int threadCount(8);
    ctpl::thread_pool pool(threadCount);

    /// Initialize all thread-local edge data:
    std::vector<EdgeThreadLocalData> edgeDataPool(threadCount);
    for (auto& edgeData : edgeDataPool)
    {
        edgeData.edgeTrackerPtr.reset(new EdgeRuntimeTracker(edgeTrackerStreamPtr));
        edgeData.svFindPtr.reset(new SVFinder(opt, readScanner, bamHeader, cset.getAllSampleReadCounts(),
            edgeData.edgeTrackerPtr, edgeData.edgeStatMan));
    }

    std::vector<std::future<GSCEdgeStats>> edgeResults;
    while (edger.next())
    {
        edgeResults.push_back(
            pool.push(
                processEdge, std::cref(opt), std::cref(cset), std::ref(edgeDataPool), edger.getEdge()));
    }

    pool.stop(true);

    GSCEdgeStats mergedStats;
    for (auto& edgeResult : edgeResults)
    {
        const GSCEdgeStats& edgeStats(edgeResult.get());
        mergedStats.merge(edgeStats);
    }

    if (! opt.edgeStatsFilename.empty())
    {
        mergedStats.save(opt.edgeStatsFilename.c_str());
    }
    if (! opt.edgeStatsReportFilename.empty())
    {
        mergedStats.report(opt.edgeStatsReportFilename.c_str());
    }
}



void
GenerateSVCandidates::
runInternal(int argc, char* argv[]) const
{
    GSCOptions opt;

    parseGSCOptions(*this,argc,argv,opt);
#ifdef DEBUG_GSV
    opt.isVerbose=true;
#endif
    runGSC(opt, name(), version());
}
