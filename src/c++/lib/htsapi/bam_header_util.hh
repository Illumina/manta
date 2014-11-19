// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Manta
// Copyright (c) 2013-2014 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/sequencing/licenses/>
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
