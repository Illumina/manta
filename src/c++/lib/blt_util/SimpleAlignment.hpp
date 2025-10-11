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

#include "blt_util/align_path.hpp"
#include "blt_util/known_pos_range2.hpp"

#include <iosfwd>

/// \brief Minimal read alignment information
///
/// Alignment information is typically processed from a BAM alignment record, but may come from other sources
struct SimpleAlignment {
  bool              is_fwd_strand = true;
  int32_t           tid           = 0;
  pos_t             pos           = 0;
  ALIGNPATH::path_t path;
};

std::ostream& operator<<(std::ostream& os, const SimpleAlignment& sa);

/// convert segment_type to match if the segment exists before or after all match segments currently in the
/// alignment
///
SimpleAlignment matchifyEdgeSegmentType(
    const SimpleAlignment&   al,
    const ALIGNPATH::align_t segmentType,
    const bool               isMatchLeadingEdge  = true,
    const bool               isMatchTrailingEdge = true);

/// transform an alignment such that any soft-clipped edge segments
/// become match.
///
/// segments are joined and start pos is adjusted appropriately
///
inline SimpleAlignment matchifyEdgeSoftClip(const SimpleAlignment& al)
{
  return matchifyEdgeSegmentType(al, ALIGNPATH::SOFT_CLIP);
}

/// get the range in reference coordinates if you did run matchifyEdgeSoftClip on an alignment:
known_pos_range2 matchifyEdgeSoftClipRefRange(const SimpleAlignment& al);
