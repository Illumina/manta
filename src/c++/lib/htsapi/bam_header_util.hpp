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

/// \brief bam record manipulation functions
///

#pragma once

#include "bam_header_info.hpp"
#include "bam_util.hpp"

#include <string>

/// parse a bam region into chrom/begin/end values
///
/// \param[out] begin_pos start position (zero-indexed, closed)
/// \param[out] end_pos end position (zero-indexed, open)
///
/// example: "chr2:100-200" will be parsed to begin_pos=99 and end_pos=200
void parse_bam_region(const char* region, std::string& chrom, int32_t& begin_pos, int32_t& end_pos);

/// parse a bam region into chrom-index/begin/end values based
/// on chromosome index lookup and end positions in bam header
///
/// \param[out] tid htslib zero-indexed contig id
/// \param[out] begin_pos start position (zero-indexed, closed)
/// \param[out] end_pos end position (zero-indexed, open)
void parse_bam_region_from_hdr(
    const bam_hdr_t* header, const char* region, int32_t& tid, int32_t& begin_pos, int32_t& end_pos);

void parse_bam_region(
    const bam_header_info& header, const char* region, int32_t& tid, int32_t& begin_pos, int32_t& end_pos);

/// return true only if the headers refer to the same
/// reference sequences in the same order.
///
bool check_header_compatibility(const bam_hdr_t& h1, const bam_hdr_t& h2);

/// try to determine the sample_name from the BAM/CRAM header
/// if none found return default string value
std::string get_bam_header_sample_name(const bam_hdr_t& header, const char* default_sample_name = "SAMPLE");
