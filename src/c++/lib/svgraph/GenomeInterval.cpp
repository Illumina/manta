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

#include "svgraph/GenomeInterval.hpp"

#include <iostream>

static std::string getChromName(const bam_header_info& bamHeader, const int tid)
{
  if (tid >= 0) {
    assert(tid < static_cast<int>(bamHeader.chrom_data.size()));
    return bamHeader.chrom_data[tid].label;
  } else {
    return "UNKNOWN";
  }
}

void summarizeGenomeInterval(const bam_header_info& bamHeader, const GenomeInterval& gi, std::ostream& os)
{
  os << "EndUserGenomeInterval: " << getChromName(bamHeader, gi.tid) << ":" << gi.range.begin_pos() + 1 << "-"
     << gi.range.end_pos();
}

std::ostream& operator<<(std::ostream& os, const GenomeInterval& gi)
{
  os << "GenomeInterval: " << gi.tid << ":" << gi.range;
  return os;
}
