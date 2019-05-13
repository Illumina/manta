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

#pragma once

#include <iosfwd>
#include <map>
#include <string>
#include <vector>

#include "blt_util/thirdparty_push.h"

#include "boost/serialization/map.hpp"
#include "boost/serialization/string.hpp"
#include "boost/serialization/vector.hpp"

#include "blt_util/thirdparty_pop.h"

#include "bam_util.hpp"

/// \brief Subset of information from a BAM file header
///
/// This holds a subset of information from the samtools bam_hdr_t object,
/// but with easier copy semantics / serialization, etc.
///
/// Currently this stores for each chromosome, the label, size and corresponding BAM index id.
///
struct bam_header_info {
  bam_header_info() = default;

  explicit bam_header_info(const bam_hdr_t& header);

  bool operator==(const bam_header_info& rhs) const
  {
    const unsigned data_size(chrom_data.size());
    if (chrom_data.size() != rhs.chrom_data.size()) return false;
    for (unsigned i(0); i < data_size; ++i) {
      if (chrom_data[i] == rhs.chrom_data[i]) continue;
      return false;
    }
    return true;
  }

  void clear()
  {
    chrom_data.clear();
    chrom_to_index.clear();
  }

  bool empty() const { return (chrom_data.empty() && chrom_to_index.empty()); }

  template <class Archive>
  void serialize(Archive& ar, const unsigned /* version */)
  {
    ar& chrom_data;
    ar& chrom_to_index;
  }

  struct chrom_info {
    explicit chrom_info(const char* init_label = nullptr, const unsigned init_length = 0)
      : label((nullptr == init_label) ? "" : init_label), length(init_length)
    {
    }

    bool operator==(const chrom_info& rhs) const { return ((label == rhs.label) && (length == rhs.length)); }

    template <class Archive>
    void serialize(Archive& ar, const unsigned /* version */)
    {
      ar& label& length;
    }

    std::string label;
    unsigned    length;
  };

  std::vector<chrom_info>        chrom_data;
  std::map<std::string, int32_t> chrom_to_index;
};

/// \brief Print header info to stream in a simple tabular format
std::ostream& operator<<(std::ostream& os, const bam_header_info& bhi);
