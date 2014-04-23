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

#pragma once

#include "blt_util/align_path.hh"

#include <iosfwd>

/// minimal summary of a query sequence aligned to a reference, roughly
/// following bam conventions for describing the alignment (apath trivially
/// maps to CIGAR string segments)
struct Alignment
{
    void
    clear()
    {
        beginPos = 0;
        apath.clear();
    }

    bool
    isAligned()
    const
    {
        return (! apath.empty());
    }

    pos_t beginPos;
    ALIGNPATH::path_t apath;
};

std::ostream&
operator<<(std::ostream& os, const Alignment& align);



struct AlignState
{
    // note the order of this enumerator is important for bit packing in client code, in particular
    // we rely on fitting the [MATCH->JUMP] states in 2 bits for the standard jump aligner
    enum index_t
    {
        MATCH,
        DELETE,
        INSERT,
        JUMP, // allows for an arbitrarily large hop between two reference regions
        SPLICE,
        JUMPINS = SPLICE, // analogous to jump state, but for very large insertions, reuse SPLICE state, in current applications we don't need both states in the same model.
        SIZE
    };

    static
    const char*
    label(const index_t i)
    {
        switch (i)
        {
        case MATCH:
            return "MATCH";
        case DELETE:
            return "DELETE";
        case INSERT:
            return "INSERT";
        case JUMP:
            return "JUMP";
        case SPLICE:
            return "SPLICE/JUMPINS";
        default:
            return "UNKNOWN";
        }
    }
};
