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

#include "depth_buffer_util.hpp"

void add_alignment_to_depth_buffer(const pos_t& pos, const ALIGNPATH::path_t& apath, depth_buffer& db)
{
  using namespace ALIGNPATH;

  pos_t ref_head_pos(pos);

  for (const path_segment& ps : apath) {
    if (is_segment_align_match(ps.type)) {
      for (unsigned j(0); j < ps.length; ++j) db.inc(ref_head_pos + static_cast<pos_t>(j));
    }

    if (is_segment_type_ref_length(ps.type)) ref_head_pos += ps.length;
  }
}
