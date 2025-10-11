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

#include "boost/any.hpp"

#include "blt_util/io_util.hpp"
#include "htsapi/bam_header_info.hpp"
#include "manta/EventInfo.hpp"
#include "manta/JunctionIdGenerator.hpp"
#include "manta/SVCandidate.hpp"
#include "manta/SVCandidateAssemblyData.hpp"
#include "manta/SVCandidateSetData.hpp"
#include "manta/SVModelScoreInfo.hpp"

#include <iosfwd>

struct VcfWriterSV {
  VcfWriterSV(
      const std::string&     referenceFilename,
      const bam_header_info& bamHeaderInfo,
      const std::string&     outputFilename,
      const bool&            isOutputContig);

  virtual ~VcfWriterSV() {}

  void writeHeader(
      const char* progName, const char* progVersion, const std::vector<std::string>& sampleNames) const;

  typedef std::vector<std::string>                                      InfoTag_t;
  typedef std::vector<std::pair<std::string, std::vector<std::string>>> SampleTag_t;

protected:
  void writeHeaderPrefix(const char* progName, const char* progVersion, std::ostream& os) const;

  void writeHeaderColumnKey(const std::vector<std::string>& sampleNames, std::ostream& os) const;

  virtual void addHeaderInfo(std::ostream& /*os*/) const {}

  virtual void addHeaderFormat(std::ostream& /*os*/) const {}

  virtual void addHeaderFilters(std::ostream& /*os*/) const {}

  void writeSVCore(
      const SVCandidateSetData&      svData,
      const SVCandidateAssemblyData& adata,
      const SVCandidate&             sv,
      const SVId&                    svId,
      const SVScoreInfo*             baseScoringInfoPtr,
      const boost::any               specializedScoringInfo,
      const EventInfo&               event,
      const bool                     isForceIntraChromBnd = false) const;

  /// add info tags which can be customized by sub-class
  virtual void modifyInfo(
      const EventInfo& /*event*/, const boost::any /*specializedScoringInfo*/, InfoTag_t& /*infotags*/) const
  {
  }

  /// add info tags specific to translocations:
  virtual void modifyTranslocInfo(
      const SVCandidate& /*sv*/,
      const SVScoreInfo* /*baseScoringInfoPtr*/,
      const bool /*isFirstOfPair*/,
      const SVCandidateAssemblyData& /*assemblyData*/,
      InfoTag_t& /*infoTags*/) const
  {
  }

  /// add info tags specific to non-translocations:
  virtual void modifyInvdelInfo(
      const SVCandidate& /*sv*/, const bool /*isBp1First*/, InfoTag_t& /*infoTags*/) const
  {
  }

  virtual void writeQual(const boost::any /*specializedScoringInfo*/, std::ostream& os) const { os << '.'; }

  virtual void writeFilter(const boost::any /*specializedScoringInfo*/, std::ostream& os) const { os << '.'; }

  virtual void modifySample(
      const SVCandidate& /*sv*/,
      const SVScoreInfo* /*baseScoringInfoPtr*/,
      const boost::any /*specializedScoringInfo*/,
      SampleTag_t& /*sampletags*/) const
  {
  }

  static void writeFilters(const std::set<std::string>& filters, std::ostream& os);

  static void writeFilters(const std::set<std::string>& filters, std::string& s);

private:
  /// \param[in] isFirstBreakend if true report bp1, else report bp2
  void writeTransloc(
      const SVCandidate&             sv,
      const SVId&                    svId,
      const SVScoreInfo*             baseScoringInfoPtr,
      const boost::any               specializedScoringInfo,
      const bool                     isFirstBreakend,
      const SVCandidateSetData&      svData,
      const SVCandidateAssemblyData& adata,
      const EventInfo&               event) const;

  void writeTranslocPair(
      const SVCandidate&             sv,
      const SVId&                    svId,
      const SVScoreInfo*             baseScoringInfoPtr,
      const boost::any               specializedScoringInfo,
      const SVCandidateSetData&      svData,
      const SVCandidateAssemblyData& adata,
      const EventInfo&               event) const;

  /// \param isIndel if true, the variant is a simple right/left breakend insert/delete combination
  void writeIndel(
      const SVCandidate& sv,
      const SVId&        svId,
      const SVScoreInfo* baseScoringInfoPtr,
      const boost::any   specializedScoringInfo,
      const bool         isIndel,
      const EventInfo&   event) const;

protected:
  const std::string&               _referenceFilename;
  const bool&                      _isOutputContig;
  mutable SynchronizedOutputStream _stream;

private:
  const bam_header_info& _header;
};
