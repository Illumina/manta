// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Manta
// Copyright (c) 2013 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/downloads/sequencing/licenses/>.
//

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
