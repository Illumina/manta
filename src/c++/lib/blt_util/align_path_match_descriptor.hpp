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
/// \brief functions to convert alignments to/from the match descriptor (MD) format
///
/// Note MD is an older format which should be considered deprecated.
///

#pragma once

#include "blt_util/align_path.hpp"

namespace ALIGNPATH {

void export_md_to_apath(
    const char* md, const bool is_fwd_strand, path_t& apath, const bool is_edge_deletion_error = true);

void apath_to_export_md(
    path_t&            apath,
    const char*        ref_seq,
    const char*        ref_end,
    const int32_t      ref_pos,
    const std::string& read_bases,
    const bool         is_fwd_strand,
    std::string&       md);

}  // namespace ALIGNPATH
