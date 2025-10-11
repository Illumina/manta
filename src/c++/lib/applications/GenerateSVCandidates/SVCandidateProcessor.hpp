//
// Manta - Structural Variant and Indel Caller
// Copyright (c) 2013-2025 Illumina, Inc.
//
// This program is licensed under the terms of the Polyform strict license
//
// ***As far as the law allows, the software comes as is, without
// any warranty or condition, and the licensor will not be liable
// to you for any damages arising out of these terms or the use
// or nature of the software, under any kind of legal claim.***
//
// You should have received a copy of the PolyForm Strict License 1.0.0
// along with this program.  If not, see <https://polyformproject.org/licenses/strict/1.0.0>.
//
//

/// \file
/// \author Chris Saunders
///

#pragma once

#include "EdgeRuntimeTracker.hpp"
#include "GSCEdgeStatsManager.hpp"
#include "GSCOptions.hpp"
#include "SVCandidateAssemblyRefiner.hpp"
#include "SVEvidenceWriter.hpp"
#include "SVScorer.hpp"
#include "SVWriter.hpp"
#include "manta/JunctionIdGenerator.hpp"
#include "manta/SVCandidateSetData.hpp"
#include "manta/SVLocusScanner.hpp"
#include "manta/SVMultiJunctionCandidate.hpp"

struct SVCandidateProcessor {
  SVCandidateProcessor(
      const GSCOptions&                           opt,
      const SVLocusScanner&                       readScanner,
      const SVLocusSet&                           cset,
      const SVWriter&                             svWriter,
      std::shared_ptr<SVEvidenceWriterSharedData> svEvidenceWriterSharedData,
      std::shared_ptr<EdgeRuntimeTracker>         edgeTrackerPtr,
      GSCEdgeStatsManager&                        edgeStatMan);

  /// Refine initial low-resolution candidates using an assembly step, then score and output final SVs
  void evaluateCandidates(
      const EdgeInfo&                              edge,
      const std::vector<SVMultiJunctionCandidate>& mjSVs,
      const SVCandidateSetData&                    svData);

private:
  void scoreSV(
      const SVCandidateSetData&                   svData,
      const std::vector<SVCandidateAssemblyData>& mjAssemblyData,
      const SVMultiJunctionCandidate&             mjSV,
      const std::vector<SVId>&                    junctionSVId,
      std::vector<bool>&                          isJunctionFiltered,
      bool&                                       isMultiJunctionEvent);

  /// Complete scoring a single SV and then write it out to VCF
  void scoreAndWriteSV(
      const EdgeInfo&                             edge,
      const SVCandidateSetData&                   svData,
      const std::vector<SVCandidateAssemblyData>& assemblyData,
      const SVMultiJunctionCandidate&             mjSV,
      const std::vector<bool>&                    isInputJunctionFiltered);

  void evaluateCandidate(
      const EdgeInfo&                 edge,
      const SVMultiJunctionCandidate& mjCandidateSV,
      const SVCandidateSetData&       svData,
      const bool                      isFindLargeInsertions);

  const GSCOptions&                   _opt;
  const SVLocusSet&                   _cset;
  const SVWriter&                     _svWriter;
  SVEvidenceWriter                    _svEvidenceWriter;
  std::shared_ptr<EdgeRuntimeTracker> _edgeTrackerPtr;
  GSCEdgeStatsManager&                _edgeStatMan;
  SVCandidateAssemblyRefiner          _svRefine;
  SVScorer                            _svScorer;
  JunctionIdGenerator                 _idgen;

  /// These are only cached here to reduce syscalls:
  std::vector<SVModelScoreInfo> _mjModelScoreInfo;
  SVModelScoreInfo              _mjJointModelScoreInfo;
  SVEvidenceWriterData          _svEvidenceWriterData;
};
