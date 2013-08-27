
///
/// \author Chris Saunders
///

#pragma once

#include "blt_util/align_path.hh"
#include "blt_util/depth_buffer.hh"


/// parse alignment into depth buffer object:
///
void
add_alignment_to_depth_buffer(
    const pos_t& pos,
    const ALIGNPATH::path_t& apath,
    depth_buffer& db);
