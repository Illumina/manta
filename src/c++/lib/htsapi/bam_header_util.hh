// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Manta - Structural Variant and Indel Caller
// Copyright (c) 2013-2015 Illumina, Inc.
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

///
/// \author Chris Saunders
///

/// \brief bam record manipulation functions
///

#pragma once

#include "bam_util.hh"
#include "bam_header_info.hh"

#include <string>


void
parse_bam_region(
    const bam_header_info& header,
    const std::string& region,
    int32_t& tid,
    int32_t& begin_pos,
    int32_t& end_pos);


/// return true only if the headers refer to the same
/// reference sequences in the same order.
///
bool
check_header_compatibility(
    const bam_header_t* h1,
    const bam_header_t* h2);


/// try to determine the sample_name from the BAM header
/// if none found return default string value
std::string
get_bam_header_sample_name(
    const std::string& bam_header_text,
    const char* default_sample_name = "SAMPLE");
