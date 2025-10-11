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
/// \author Xiaoyu Chen
///

#pragma once

#include "format/VcfWriterSV.hpp"
#include "options/CallOptionsTumor.hpp"

struct VcfWriterTumorSV : public VcfWriterSV {
  VcfWriterTumorSV(
      const CallOptionsTumor& tumorOpt,
      const bool              isMaxDepthFilter,
      const std::string&      referenceFilename,
      const bam_header_info&  bamHeaderInfo,
      const std::string&      outputFilename,
      const bool&             isOutputContig)
    : VcfWriterSV(referenceFilename, bamHeaderInfo, outputFilename, isOutputContig),
      _tumorOpt(tumorOpt),
      _isMaxDepthFilter(isMaxDepthFilter)
  {
  }

  void writeSV(
      const SVCandidateSetData&      svData,
      const SVCandidateAssemblyData& adata,
      const SVCandidate&             sv,
      const SVId&                    svId,
      const SVScoreInfo&             baseScoringInfo,
      const SVScoreInfoTumor&        tumorInfo,
      const EventInfo&               event) const;

private:
  void addHeaderInfo(std::ostream& os) const override;

  void addHeaderFormat(std::ostream& os) const override;

  void addHeaderFilters(std::ostream& os) const override;

  void writeFilter(const boost::any specializedScoringInfo, std::ostream& os) const override;

  void modifySample(
      const SVCandidate& sv,
      const SVScoreInfo* baseScoringInfoPtr,
      const boost::any   specializedScoringInfo,
      SampleTag_t&       sampletags) const override;

  void modifyTranslocInfo(
      const SVCandidate&             sv,
      const SVScoreInfo*             baseScoringInfoPtr,
      const bool                     isFirstOfPair,
      const SVCandidateAssemblyData& assemblyData,
      InfoTag_t&                     infotags) const override;

  const CallOptionsTumor& _tumorOpt;
  const bool              _isMaxDepthFilter;
};
