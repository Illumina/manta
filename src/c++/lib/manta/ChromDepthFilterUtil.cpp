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

#include "manta/ChromDepthFilterUtil.hpp"
#include "common/Exceptions.hpp"

#include <sstream>

ChromDepthFilterUtil::ChromDepthFilterUtil(
    const std::string& chromDepthFile, const double maxDepthFactor, const bam_header_info& header)
  : _isMaxDepthFilter(!chromDepthFile.empty())
{
  using namespace illumina::common;

  // read in chrom depth file if one is specified:
  if (!_isMaxDepthFilter) return;

  cdmap_t chromDepth;
  parse_chrom_depth(chromDepthFile, chromDepth);

  // translate string chrom labels into tid values in lookup vector:
  //
  unsigned int callableChromCount(0);
  for (const bam_header_info::chrom_info& cdata : header.chrom_data) {
    cdmap_t::const_iterator cdi(chromDepth.find(cdata.label));

    if (cdi != chromDepth.end()) {
      _maxDepthFilter.push_back(cdi->second * maxDepthFactor);
      callableChromCount++;
    } else {
      _maxDepthFilter.push_back(0);
    }
    assert(_maxDepthFilter.back() >= 0.);
  }

  if (callableChromCount != chromDepth.size()) {
    std::ostringstream oss;
    oss << chromDepth.size() << " chromosomes in chrom depth file, but " << callableChromCount
        << " found in the bam header";
    BOOST_THROW_EXCEPTION(GeneralException(oss.str()));
  }
}
