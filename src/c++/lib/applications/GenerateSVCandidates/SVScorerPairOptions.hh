// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Manta
// Copyright (c) 2013-2015 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/sequencing/licenses/>
//

///
/// \author Chris Saunders and Xiaoyu Chen
///

#pragma once

#include "blt_util/blt_types.hh"


/// shared options related to read pair support:
struct PairOptions
{
    PairOptions(const bool isRNA) :
        minFragSupport(50),
        minFragProb(! isRNA ? 0.0001 : 0.0)
    {}

    /// we're interested in any fragments which cross center pos with at least N bases of support on each side
    /// (note this definition is certain to overlap the split read definition whenever N is less than the read length
    ///
    /// for reads shorter than this length, the whole read is required...
    const pos_t minFragSupport;

    const float minFragProb;
};
