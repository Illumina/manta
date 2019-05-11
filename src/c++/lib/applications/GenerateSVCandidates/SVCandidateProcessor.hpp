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
