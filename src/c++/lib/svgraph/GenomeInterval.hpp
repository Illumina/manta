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

#include "blt_util/known_pos_range2.hpp"
#include "htsapi/bam_header_info.hpp"

///
/// \author Chris Saunders
///
/// \brief GenomeInterval identifies a contiguous chromosomal region.
///
/// \details GenomeInterval identifies single contiguous chromosome range. All internal locations use a
/// chromosome index number. GenomeInterval uses boost::serialize to save/load the class.
struct GenomeInterval {
  GenomeInterval(const int32_t initTid = 0, const pos_t beginPos = 0, const pos_t endPos = 0)
    : tid(initTid), range(beginPos, endPos)
  {
  }

  /// \brief Identify if the GenomeIntersect Intersects with another GenomeInterval
  ///
  /// 1. The ids must be the same
  /// 2. The range of the GenomeIntervals must overlap
  bool isIntersect(const GenomeInterval& gi) const
  {
    if (tid != gi.tid) return false;
    return range.is_range_intersect(gi.range);
  }

  bool operator<(const GenomeInterval& rhs) const
  {
    if (tid < rhs.tid) return true;
    if (tid == rhs.tid) {
      return (range < rhs.range);
    }
    return false;
  }

  bool operator==(const GenomeInterval& rhs) const { return ((tid == rhs.tid) && (range == rhs.range)); }

  void clear()
  {
    tid = 0;
    range.clear();
  }

  template <class Archive>
  void serialize(Archive& ar, const unsigned /* version */)
  {
    ar& tid& range;
  }

  /// \brief Samtools chromosome index (samtools)
  int32_t tid;
  /// \brief Chromosome interval range
  known_pos_range2 range;
};

/// Pretty print summary information from a genome interval for end-user error message
void summarizeGenomeInterval(const bam_header_info& bamHeader, const GenomeInterval& gi, std::ostream& os);

/// Debug printer for genome interval
std::ostream& operator<<(std::ostream& os, const GenomeInterval& gi);

BOOST_CLASS_IMPLEMENTATION(GenomeInterval, boost::serialization::object_serializable)
