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

#include "bam_header_info.hpp"

#include <iostream>

bam_header_info::bam_header_info(const bam_hdr_t& header)
{
  for (int i(0); i < header.n_targets; ++i) {
    chrom_data.emplace_back(header.target_name[i], header.target_len[i]);
    chrom_to_index[header.target_name[i]] = (int32_t)i;
  }
}

std::ostream& operator<<(std::ostream& os, const bam_header_info& bhi)
{
  unsigned index(0);

  os << "chomosome_id_map:\n";
  for (const bam_header_info::chrom_info& info : bhi.chrom_data) {
    os << "index: " << index << " label: " << info.label << " length: " << info.length << '\n';
    index++;
  }
  return os;
}
