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

/// \author Chris Saunders
///

#pragma once

#include "blt_util/bam_util.hh"

#include <string>
#include <vector>


/// minimal c++ bam header info
///
/// this class replicates the minimum information
/// from the bam header required to parse regions
/// (ie. chr1:20-30). It is friendlier to mem management
/// and serialization than using the samtools struct
///
struct bam_header_info
{

    bam_header_info()
    {}

    bam_header_info(const bam_header_t& header);

    template<class Archive>
    void serialize(Archive & ar, const unsigned /* version */)
    {
        ar & chrom_data;
    }

    struct chrom_info
    {
        chrom_info(
                const char* init_label,
                const unsigned init_length) :
                    label(init_label),
                    length(init_length)
        {}

        std::string label;
        unsigned length;
    };

    std::vector<chrom_info> chrom_data;
};
