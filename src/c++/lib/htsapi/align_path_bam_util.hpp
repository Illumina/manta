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
#include "htsapi/bam_util.hpp"

/// Convert htslib CIGAR format into ALIGNPATH::path_t
///
void bam_cigar_to_apath(const uint32_t* bam_cigar, const unsigned n_cigar, ALIGNPATH::path_t& apath);

/// Convert ALIGNPATH::path_t into htslib CIGAR format
///
/// bam_cigar must already be allocated to at least apath.size() capacity
///
void apath_to_bam_cigar(const ALIGNPATH::path_t& apath, uint32_t* bam_cigar);

/// Convert ALIGNPATH::path_t into htslib CIGAR format, and insert this
/// into the input htslib BAM record
///
void edit_bam_cigar(const ALIGNPATH::path_t& apath, bam1_t& br);
