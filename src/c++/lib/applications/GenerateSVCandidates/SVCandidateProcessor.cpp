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
/// \author Naoki Nariai
///

#include "SVCandidateProcessor.hpp"

#include <iostream>

#include "blt_util/log.hpp"
#include "manta/SVMultiJunctionCandidateUtil.hpp"
#include "svgraph/EdgeInfoUtil.hpp"
#include "svgraph/SVLocusSet.hpp"

SVCandidateProcessor::SVCandidateProcessor(
    const GSCOptions&                           opt,
    const SVLocusScanner&                       readScanner,
    const SVLocusSet&                           cset,
    const SVWriter&                             svWriter,
    std::shared_ptr<SVEvidenceWriterSharedData> svEvidenceWriterSharedData,
    std::shared_ptr<EdgeRuntimeTracker>         edgeTrackerPtr,
    GSCEdgeStatsManager&                        edgeStatMan)
  : _opt(opt),
    _cset(cset),
    _svWriter(svWriter),
    _svEvidenceWriter(opt, svEvidenceWriterSharedData),
    _edgeTrackerPtr(edgeTrackerPtr),
    _edgeStatMan(edgeStatMan),
    _svRefine(opt, cset.getBamHeader(), cset.getAllSampleReadCounts(), _edgeTrackerPtr),
    _svScorer(opt, readScanner, cset.getBamHeader()),
    _svEvidenceWriterData(opt.alignFileOpt.alignmentFilenames.size())
{
}

void SVCandidateProcessor::evaluateCandidates(
    const EdgeInfo&                              edge,
    const std::vector<SVMultiJunctionCandidate>& mjSVs,
    const SVCandidateSetData&                    svData)
{
  const bool isIsolatedEdge(testIsolatedEdge(_cset, edge));

  bool isFindLargeInsertions(isIsolatedEdge);
  if (isFindLargeInsertions) {
    for (const SVMultiJunctionCandidate& mjCandidateSV : mjSVs) {
      for (const SVCandidate& candidateSV : mjCandidateSV.junction) {
        if (!isComplexSV(candidateSV)) isFindLargeInsertions = false;
      }
    }
  }

  _svRefine.clearEdgeData();
  _svEvidenceWriterData.clear();

  for (const auto& cand : mjSVs) {
    evaluateCandidate(edge, cand, svData, isFindLargeInsertions);
  }

  _svEvidenceWriter.write(_svEvidenceWriterData);
}

static bool isAnyFalse(const std::vector<bool>& vb)
{
  for (const bool val : vb) {
    if (!val) return true;
  }
  return false;
}

/// Each junction of a candidate SV is checked and marked as filtered if it is found to be low quality.
///
/// The junction filtration rules are:
/// 1. If the junction is part of an sv candidate where all junctions have a spanning evidence count less than
/// minCandidateSpanningCount, then all junctions in the sv candidate will be filtered.
/// 2. Any individual junction with a spanning evidence count less than minJunctionCandidateSpanningCount will
/// be filtered.
/// 3. Complex candidate regions which did not produce a successful contig alignment will be filtered (the
/// candidate is not spanning so there is no imprecise hypothesis to pursue if contig alignment fails).
/// 4. Candidates will be filtered if their size is smaller than minCandidateVariantSize
///
static void checkJunctionsToFilter(
    const SVMultiJunctionCandidate&             mjSV,
    const std::vector<SVCandidateAssemblyData>& mjAssemblyData,
    std::vector<bool>&                          isJunctionFiltered,
    const GSCOptions&                           opt)
{
  const unsigned junctionCount(mjSV.junction.size());

  // relax min spanning count to just 2 for the special case of a
  // junction that is part of a multi-junction EVENT
  const unsigned minJunctionCandidateSpanningCount(std::min(2u, opt.minCandidateSpanningCount));

  // first step of SV filtering (before candidates are written)
  // 2 junction filter types:
  // 1) tests where the junction can fail independently
  // 2) tests where all junctions have to fail for the candidate to be filtered:
  bool isCandidateSpanFail(true);

  for (unsigned junctionIndex(0); junctionIndex < junctionCount; ++junctionIndex) {
    const SVCandidateAssemblyData& assemblyData(mjAssemblyData[junctionIndex]);
    const SVCandidate&             sv(mjSV.junction[junctionIndex]);

    const bool isCandidateSpanning(assemblyData.isCandidateSpanning);

#ifdef DEBUG_GSV
    log_os << __FUNCTION__ << ": isSpanningSV junction: " << isCandidateSpanning << "\n";
    log_os << __FUNCTION__ << ": Pairs:" << sv.bp1.getPairCount() << " Spanning:" << sv.bp1.getSpanningCount()
           << "\n";
#endif

    // junction dependent tests:
    //   (1) at least one junction in the set must have spanning count of minCandidateSpanningCount or more
    bool isJunctionSpanFail(false);
    if (isCandidateSpanning) {
      if (sv.getPostAssemblySpanningCount(opt.isRNA) < opt.minCandidateSpanningCount) {
        isJunctionSpanFail = true;
      }
    }
    if (!isJunctionSpanFail) isCandidateSpanFail = false;

    // independent tests -- as soon as one of these fails, we can continue:
    if (isCandidateSpanning) {
      //   (2) each spanning junction in the set must have spanning count of
      //      minJunctionCandidateSpanningCount or more
      if (sv.getPostAssemblySpanningCount(opt.isRNA) < minJunctionCandidateSpanningCount) {
        isJunctionFiltered[junctionIndex] = true;
        continue;
      }
    } else {
      if (sv.isImprecise()) {
        // (3) no unassembled non-spanning candidates, in this case a non-spanning low-res candidate
        // went into assembly but did not produce a successful contig alignment:
#ifdef DEBUG_GSV
        log_os << __FUNCTION__ << ": Rejecting candidate junction: imprecise non-spanning SV\n";
#endif
        isJunctionFiltered[junctionIndex] = true;
        continue;
      }
    }

    // (4) check min size, or if it is IMPRECISE case where CIEND is a subset of CIPOS for candidate output:
    if (isSVBelowMinSize(sv, opt.scanOpt.minCandidateVariantSize)) {
#ifdef DEBUG_GSV
      log_os << __FUNCTION__ << ": Filtering out candidate below min size before candidate output stage\n";
#endif
      isJunctionFiltered[junctionIndex] = true;
      continue;
    }
  }

  if (isCandidateSpanFail) {
    for (unsigned junctionIndex(0); junctionIndex < junctionCount; ++junctionIndex) {
#ifdef DEBUG_GSV
      log_os << __FUNCTION__ << ": Rejecting candidate junction: minCandidateSpanningCount\n";
#endif
      isJunctionFiltered[junctionIndex] = true;
    }
  }
}

void SVCandidateProcessor::scoreSV(
    const SVCandidateSetData&                   svData,
    const std::vector<SVCandidateAssemblyData>& mjAssemblyData,
    const SVMultiJunctionCandidate&             mjSV,
    const std::vector<SVId>&                    junctionSVId,
    std::vector<bool>&                          isJunctionFiltered,
    bool&                                       isMJEvent)
{
  if (_opt.isSkipScoring) return;

  // check min size for scoring:
  const unsigned junctionCount(isJunctionFiltered.size());
  for (unsigned junctionIndex(0); junctionIndex < junctionCount; ++junctionIndex) {
    if (isJunctionFiltered[junctionIndex]) continue;

    const SVCandidate& sv(mjSV.junction[junctionIndex]);
    if (isSVBelowMinSize(sv, _opt.minScoredVariantSize)) {
#ifdef DEBUG_GSV
      log_os << __FUNCTION__ << ": Filtering out candidate junction below min size at scoring stage\n";
#endif
      isJunctionFiltered[junctionIndex] = true;
    }
  }

  // check to see if all junctions are filtered before scoring:
  //
  if (!isAnyFalse(isJunctionFiltered)) return;

  _svScorer.scoreSV(
      svData,
      mjAssemblyData,
      mjSV,
      junctionSVId,
      isJunctionFiltered,
      _opt.isSomatic(),
      _opt.isTumorOnly(),
      _mjModelScoreInfo,
      _mjJointModelScoreInfo,
      isMJEvent,
      _svEvidenceWriterData);
}

void SVCandidateProcessor::scoreAndWriteSV(
    const EdgeInfo&                             edge,
    const SVCandidateSetData&                   svData,
    const std::vector<SVCandidateAssemblyData>& mjAssemblyData,
    const SVMultiJunctionCandidate&             mjSV,
    const std::vector<bool>&                    isInputJunctionFiltered)
{
  // track filtration for each junction:
  std::vector<bool> isCandidateJunctionFiltered(isInputJunctionFiltered);
  checkJunctionsToFilter(mjSV, mjAssemblyData, isCandidateJunctionFiltered, _opt);

  // check to see if all junctions are filtered, if so skip the whole candidate:
  //
  if (!isAnyFalse(isCandidateJunctionFiltered)) {
#ifdef DEBUG_GSV
    log_os << __FUNCTION__ << ": Rejecting candidate, all junctions filtered.\n";
#endif
    return;
  }

  // compute all junction ids
  //
  const unsigned    junctionCount(isCandidateJunctionFiltered.size());
  std::vector<SVId> junctionSVId(junctionCount);
  for (unsigned junctionIndex(0); junctionIndex < junctionCount; ++junctionIndex) {
    const SVCandidate& sv(mjSV.junction[junctionIndex]);
    SVId&              svId(junctionSVId[junctionIndex]);
    _idgen.getId(edge, sv, _opt.isRNA, svId);
  }

  // Now score the SV
  //
  bool              isMultiJunctionEvent(false);
  std::vector<bool> isScoredJunctionFiltered(isCandidateJunctionFiltered);
  scoreSV(svData, mjAssemblyData, mjSV, junctionSVId, isScoredJunctionFiltered, isMultiJunctionEvent);

  // finally, write out to all VCF streams
  //
  _svWriter.writeSV(
      svData,
      mjAssemblyData,
      mjSV,
      isCandidateJunctionFiltered,
      isScoredJunctionFiltered,
      junctionSVId,
      _mjModelScoreInfo,
      _mjJointModelScoreInfo,
      isMultiJunctionEvent);
}

void SVCandidateProcessor::evaluateCandidate(
    const EdgeInfo&                 edge,
    const SVMultiJunctionCandidate& mjCandidateSV,
    const SVCandidateSetData&       svData,
    const bool                      isFindLargeInsertions)
{
  assert(!mjCandidateSV.junction.empty());

  const unsigned junctionCount(mjCandidateSV.junction.size());

  if (_opt.isVerbose) {
    log_os << __FUNCTION__ << ": Starting analysis for SV candidate containing " << junctionCount
           << " junctions. Low-resolution junction candidate ids:";
    for (const SVCandidate& sv : mjCandidateSV.junction) {
      log_os << " " << sv.candidateIndex;
    }
    log_os << "\n";
  }
#ifdef DEBUG_GSV
  log_os << __FUNCTION__ << ": CandidateSV: " << mjCandidateSV << "\n";
#endif

  const bool isComplex(isComplexSV(mjCandidateSV));
  _edgeTrackerPtr->addCandidate(isComplex);

  _edgeStatMan.updateJunctionCandidateCounts(edge, junctionCount, isComplex);

  // true if any junction returns a complex (indel) assembly result
  //
  // if true, and if the candidate is multi-junction then score and write each junction independently:
  bool                                 isAnySmallAssembler(false);
  std::vector<SVCandidateAssemblyData> mjAssemblyData(junctionCount);

  if (!_opt.isSkipAssembly) {
    const TimeScoper assmTime(_edgeTrackerPtr->assemblyTime);
    for (unsigned junctionIndex(0); junctionIndex < junctionCount; ++junctionIndex) {
      const SVCandidate&       candidateSV(mjCandidateSV.junction[junctionIndex]);
      SVCandidateAssemblyData& assemblyData(mjAssemblyData[junctionIndex]);
      try {
        _svRefine.getCandidateAssemblyData(candidateSV, isFindLargeInsertions, assemblyData);
      } catch (illumina::common::ExceptionData& e) {
        std::ostringstream oss;
        oss << "Exception caught while attempting to assemble " << candidateSV;
        e << boost::error_info<struct assembly_candidate_info, std::string>(oss.str());
        throw;
      } catch (...) {
        log_os << "Exception caught while attempting to assemble " << candidateSV << "\n";
        throw;
      }

      if (_opt.isVerbose) {
        log_os << __FUNCTION__ << ": Candidate assembly complete for junction " << junctionIndex << "/"
               << junctionCount << ". Assembled candidate count: " << assemblyData.svs.size() << "\n";
      }

      if (!assemblyData.svs.empty()) {
        const unsigned assemblyCount(assemblyData.svs.size());

        const bool isSpanning(assemblyData.isSpanning);

        _edgeStatMan.updateAssemblyCount(edge, assemblyCount, isSpanning);

        if (isSpanning) {
          // can't be multi-junction and multi-assembly at the same time:
          assert(!((junctionCount > 1) && (assemblyCount > 1)));
        } else {
          isAnySmallAssembler = true;
        }

        // fill in assembly tracking data:
        _edgeTrackerPtr->addAssembledCandidate(isComplex);
      } else {
        _edgeStatMan.updateAssemblyCount(edge, 0, assemblyData.isSpanning, assemblyData.isOverlapSkip);
      }
    }
  }

  SVMultiJunctionCandidate mjAssembledCandidateSV;
  mjAssembledCandidateSV.junction.resize(junctionCount);
  std::vector<bool> isJunctionFiltered(junctionCount, false);

  std::vector<unsigned> junctionTracker(junctionCount, 0);
  while (true) {
    // Note this loop is an accident -- it was intended to enumerate all assembly
    // combinations for  multiple junctions with multiple assemblies each.
    // It doesn't do that -- but the broken thing it does, in fact, do, is what we
    // want for the isAnySmallAssembler case so it's well enough for now.
    //
    /// TODO: update docs -- what is it we "want for the isAnySmallAssembler case"?
    bool isWrite(false);
    for (unsigned junctionIndex(0); junctionIndex < junctionCount; ++junctionIndex) {
      const SVCandidateAssemblyData& assemblyData(mjAssemblyData[junctionIndex]);
      unsigned&                      assemblyIndex(junctionTracker[junctionIndex]);

      if (assemblyData.svs.empty()) {
        if (assemblyIndex != 0) continue;
#ifdef DEBUG_GSV
        log_os << __FUNCTION__ << ": score and output low-res candidate junction " << junctionIndex << "\n";
#endif
        mjAssembledCandidateSV.junction[junctionIndex] = mjCandidateSV.junction[junctionIndex];
      } else {
        if (assemblyIndex >= assemblyData.svs.size()) continue;
        const SVCandidate& assembledSV(assemblyData.svs[assemblyIndex]);
#ifdef DEBUG_GSV
        log_os << __FUNCTION__ << ": score and output assembly candidate junction " << junctionIndex << ": "
               << assembledSV << "\n";
#endif
        mjAssembledCandidateSV.junction[junctionIndex] = assembledSV;
      }
      assemblyIndex++;
      isWrite = true;
    }
    if (!isWrite) break;

    // if any junctions go into the small assembler (for instance b/c the breakends are too close), then
    // treat all junctions as independent:
    //
    const TimeScoper scoreTime(_edgeTrackerPtr->scoreTime);
    if ((junctionCount > 1) && isAnySmallAssembler) {
      // call each junction independently by providing a filter vector in each iteration
      // which only leaves a single junction unfiltered:
      for (unsigned junctionIndex(0); junctionIndex < junctionCount; ++junctionIndex) {
        std::vector<bool> isJunctionFilteredHack(junctionCount, true);
        isJunctionFilteredHack[junctionIndex] = false;
        scoreAndWriteSV(edge, svData, mjAssemblyData, mjAssembledCandidateSV, isJunctionFilteredHack);
      }
    } else {
      scoreAndWriteSV(edge, svData, mjAssemblyData, mjAssembledCandidateSV, isJunctionFiltered);
    }
  }
}
