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
// <https://github.com/sequencing/licenses/>
//

/// \author Chris Saunders
///

#include "blt_util/bam_header_info.hh"

#include "boost/foreach.hpp"

#include <iostream>



bam_header_info::
bam_header_info(const bam_header_t& header)
{
    for (int i(0); i<header.n_targets; ++i)
    {
        chrom_data.push_back(chrom_info(header.target_name[i],header.target_len[i]));
        chrom_to_int[header.target_name[i]] = (int32_t)i;
    }
}


std::ostream&
operator<<(std::ostream& os, const bam_header_info& bhi)
{
    unsigned index(0);

    os << "chomosome_id_map:\n";
    BOOST_FOREACH(const bam_header_info::chrom_info& info, bhi.chrom_data)
    {
        os << "index: " << index << " label: " << info.label << " length: " << info.length << '\n';
        index++;
    }
    return os;
}
