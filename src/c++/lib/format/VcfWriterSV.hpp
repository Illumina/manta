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
