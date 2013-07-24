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
/// \brief align a contig across two breakend regions
///


#pragma once

#include "blt_util/align_path.hh"


struct AlignerUtil
{
    static
    void
    updatePath(
        ALIGNPATH::path_t& path,
        ALIGNPATH::path_segment& ps,
        ALIGNPATH::align_t atype)
    {
        if (ps.type == atype) return;
        if (ps.type != ALIGNPATH::NONE) path.push_back(ps);
        ps.type = atype;
        ps.length = 0;
    }
};

