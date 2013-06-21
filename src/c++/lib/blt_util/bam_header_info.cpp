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

#include "blt_util/bam_header_info.hh"


bam_header_info::
bam_header_info(const bam_header_t& header)
{
    for(int i(0);i<header.n_targets;++i)
    {
        chrom_data.push_back(chrom_info(header.target_name[i],header.target_len[i]));
    }
}
