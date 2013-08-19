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
/// \brief functions to convert alignments to/from the match descriptor (MD) format
///
/// Note MD is an older format which should be considered depricated.
///

#pragma once

#include "blt_util/align_path.hh"



namespace ALIGNPATH
{

void
export_md_to_apath(const char* md,
                   const bool is_fwd_strand,
                   path_t& apath,
                   const bool is_edge_deletion_error=true);

void
apath_to_export_md(path_t& apath,
                   const char* ref_seq,
                   const char* ref_end,
                   const int32_t ref_pos,
                   const std::string& read_bases,
                   const bool is_fwd_strand,
                   std::string& md);

}
