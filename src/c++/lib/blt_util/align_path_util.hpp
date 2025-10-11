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

#pragma once

#include "blt_util/align_path.hpp"
#include "blt_util/blt_types.hpp"

#include <cassert>

namespace ALIGNPATH {

inline void increment_path(const path_t& path, unsigned& path_index, unsigned& read_offset, pos_t& ref_offset)
{
  const path_segment& ps(path[path_index]);

  if (is_segment_align_match(ps.type)) {
    read_offset += ps.length;
    ref_offset += ps.length;
  } else if ((ps.type == DELETE) || (ps.type == SKIP)) {
    ref_offset += ps.length;
  } else if ((ps.type == INSERT) || (ps.type == SOFT_CLIP)) {
    read_offset += ps.length;
  } else if ((ps.type == HARD_CLIP) || (ps.type == PAD)) {
    // do nothing
  } else {
    assert(false && "Unexpected alignment type");  // can't handle other CIGAR types yet
  }

  path_index++;
}

// Initialize to the segment count, insert and delete size of a
// swap in the path. assumes path_index points to the begining of
// a swap:
//
struct swap_info {
  swap_info(const path_t& path, const unsigned path_index)
    : n_seg(path_index), insert_length(0), delete_length(0)
  {
    const unsigned aps(path.size());
    for (; (n_seg < aps) && is_segment_type_indel(path[n_seg].type); ++n_seg) {
      const path_segment& ps(path[n_seg]);
      if (ps.type == INSERT) {
        insert_length += ps.length;
      } else if (ps.type == DELETE) {
        delete_length += ps.length;
      } else {
        assert(false && "Unexpected alignment type");
      }
    }
    n_seg -= path_index;
  }

  unsigned n_seg;
  unsigned insert_length;
  unsigned delete_length;
};
}  // namespace ALIGNPATH
