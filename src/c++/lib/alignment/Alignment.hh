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
        alignStart=0;
        apath.clear();
    }

    unsigned alignStart;
    ALIGNPATH::path_t apath;
};
