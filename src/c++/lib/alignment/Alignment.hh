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

/// derived from ELAND implementation by Tony Cox

#pragma once

#include "blt_util/align_path.hh"


/// minimal summary of a query sequence aligned to a reference, roughly
/// following bam conventions
struct Alignment
{
    void
    clear()
    {
        alignStart = 0;
        alignEnd   = 0;
        apath.clear();
    }

    bool
    isAligned()
    const
    {
        return !apath.empty();
    }

    unsigned alignStart;
    unsigned alignEnd;
    ALIGNPATH::path_t apath;
};

struct AlignState
{
    enum index_t
    {
        MATCH,
        DELETE,
        INSERT,
        JUMP,
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
        default:
            return "UNKNOWN";
        }
    }
};
