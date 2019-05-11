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
