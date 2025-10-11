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
