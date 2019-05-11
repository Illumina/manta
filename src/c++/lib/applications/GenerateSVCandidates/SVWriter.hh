//
// Manta - Structural Variant and Indel Caller
// Copyright (c) 2013-2018 Illumina, Inc.
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

#include "GSCOptions.hh"

#include "common/OutStream.hh"
#include "format/VcfWriterCandidateSV.hh"
#include "format/VcfWriterDiploidSV.hh"
#include "format/VcfWriterRnaSV.hh"
#include "format/VcfWriterSomaticSV.hh"
#include "format/VcfWriterTumorSV.hh"
#include "manta/JunctionIdGenerator.hh"
#include "manta/SVCandidateAssemblyData.hh"
#include "manta/SVCandidateSetData.hh"
#include "manta/SVMultiJunctionCandidate.hh"

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
