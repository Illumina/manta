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
/// \brief Utility functions on BAM records

#pragma once

#include "blt_util/SimpleAlignment.hpp"
#include "blt_util/align_path.hpp"
#include "htsapi/bam_header_info.hpp"
#include "htsapi/bam_record.hpp"

#include <iosfwd>

/// \brief Test if this read is part of a pair where both members are mapped.
bool is_mapped_pair(const bam_record& bam_read);

/// \brief Test if this read part of a pair where both members are mapped to the same chrom.
///
bool is_mapped_chrom_pair(const bam_record& bam_read);

/// \brief Test if this read is part of mapped pair with "innie" orientation.
///
/// This does not test MAPQ or fragment size, but could be used as the core of a 'proper-pair' predicate.
///
/// This is designed to return true for the common case f pos == mate_pos occurring for short FFPE fragments.
bool is_innie_pair(const bam_record& bam_read);

/// Based on pair alignment information (MC tag) if available, detect cases where a read extends
/// past the beginning of its mate, suggesting the presence of unfiltered adapter sequence
/// If MC information is not available, use a more aggressive heuristic (any 3' softclipping)
/// but still less aggressive then is_possible_adapter_pair
/// Used for read filtration of assembly to avoid including adapter sequence, but also missing any evidence
/// reads Assumes innie pairs
bool is_adapter_pair(const bam_record& bamRead);

/// \brief Detect cases where paired-end reads overlap
///
/// Overlapping read pairs implies that the sequenced DNA/RNA fragment size is shorter than
/// the combined read-length,
///
/// This is an approximation because it's based on a single bam record, and as such
/// the length and alignment of the mate read are inferred. An exact answer would
/// require having BAM records for both ends of the read pair.
bool is_overlapping_pair(const bam_record& bam_read, const SimpleAlignment& matchedAlignment);

/// \return The average basecall quality score for this read.
unsigned get_avg_quality(const bam_record& bam_read);

/// select 'first' read in pair such that you
/// consistently get only one read per-pair
/// (assuming the bam file is properly formatted)
inline bool isFirstRead(const bam_record& bamRead)
{
  if (bamRead.pos() < bamRead.mate_pos()) return true;
  if ((bamRead.pos() == bamRead.mate_pos()) && bamRead.is_first()) return true;
  return false;
}

/// \return The BAM record's auxillary RG tag value, or an empty string if no RG tag exists.
inline const char* getReadGroup(const bam_record& bamRead)
{
  static const char defaultRG[] = "";
  static const char rgTag[]     = {'R', 'G'};

  const char* rgStr(bamRead.get_string_tag(rgTag));

  return ((nullptr == rgStr) ? defaultRG : rgStr);
}

/// Pretty print summary information from an alignment record for an end-user error message
void summarizeAlignmentRecord(const bam_header_info& bamHeader, const bam_record& bamRead, std::ostream& os);
