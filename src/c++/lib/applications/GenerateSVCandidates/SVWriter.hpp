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

#include "GSCOptions.hpp"
#include "common/OutStream.hpp"
#include "format/VcfWriterCandidateSV.hpp"
#include "format/VcfWriterDiploidSV.hpp"
#include "format/VcfWriterRnaSV.hpp"
#include "format/VcfWriterSomaticSV.hpp"
#include "format/VcfWriterTumorSV.hpp"
#include "manta/JunctionIdGenerator.hpp"
#include "manta/SVCandidateAssemblyData.hpp"
#include "manta/SVCandidateSetData.hpp"
#include "manta/SVMultiJunctionCandidate.hpp"

struct SVWriter {
  SVWriter(
      const GSCOptions&      initOpt,
      const bam_header_info& bamHeaderInfo,
      const char*            progName,
      const char*            progVersion);

  void writeSV(
      const SVCandidateSetData&                   svData,
      const std::vector<SVCandidateAssemblyData>& mjAssemblyData,
      const SVMultiJunctionCandidate&             mjSV,
      const std::vector<bool>&                    isCandidateJunctionFiltered,
      const std::vector<bool>&                    isScoredJunctionFiltered,
      const std::vector<SVId>&                    junctionSVId,
      const std::vector<SVModelScoreInfo>&        mjModelScoreInfo,
      const SVModelScoreInfo&                     mjJointModelScoreInfo,
      const bool                                  isMJEvent) const;

private:
  ///////////////////////// data:
  const GSCOptions& opt;
  unsigned          diploidSampleCount;

  VcfWriterCandidateSV                candWriter;
  std::unique_ptr<VcfWriterDiploidSV> diploidWriter;
  std::unique_ptr<VcfWriterSomaticSV> somWriter;
  std::unique_ptr<VcfWriterTumorSV>   tumorWriter;
  std::unique_ptr<VcfWriterRnaSV>     rnaWriter;
};
