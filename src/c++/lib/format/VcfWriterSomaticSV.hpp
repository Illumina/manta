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

#include "format/VcfWriterSV.hpp"
#include "manta/SVModelScoreInfo.hpp"
#include "options/CallOptionsSomatic.hpp"

struct VcfWriterSomaticSV : public VcfWriterSV {
  VcfWriterSomaticSV(
      const CallOptionsSomatic& somaticOpt,
      const bool                isMaxDepthFilter,
      const std::string&        referenceFilename,
      const bam_header_info&    bamHeaderInfo,
      const std::string&        outputFilename,
      const bool&               isOutputContig)
    : VcfWriterSV(referenceFilename, bamHeaderInfo, outputFilename, isOutputContig),
      _somaticOpt(somaticOpt),
      _isMaxDepthFilter(isMaxDepthFilter)
  {
  }

  void writeSV(
      const SVCandidateSetData&      svData,
      const SVCandidateAssemblyData& adata,
      const SVCandidate&             sv,
      const SVId&                    svId,
      const SVScoreInfo&             baseScoringInfo,
      const SVScoreInfoSomatic&      somaticInfo,
      const EventInfo&               event,
      const SVScoreInfoSomatic&      singleJunctionSomaticInfo) const;

private:
  void addHeaderInfo(std::ostream& os) const override;

  void addHeaderFormat(std::ostream& os) const override;

  void addHeaderFilters(std::ostream& os) const override;

  void modifyInfo(
      const EventInfo&          event,
      const boost::any          specializedScoringInfo,
      std::vector<std::string>& infotags) const override;

  void modifyTranslocInfo(
      const SVCandidate&             sv,
      const SVScoreInfo*             baseScoringInfoPtr,
      const bool                     isFirstOfPair,
      const SVCandidateAssemblyData& assemblyData,
      std::vector<std::string>&      infotags) const override;

  void modifySample(
      const SVCandidate& sv,
      const SVScoreInfo* baseScoringInfoPtr,
      const boost::any   specializedScoringInfo,
      SampleTag_t&       sampletags) const override;

  void writeFilter(const boost::any specializedScoringInfo, std::ostream& os) const override;

  const CallOptionsSomatic& _somaticOpt;
  const bool                _isMaxDepthFilter;
};
