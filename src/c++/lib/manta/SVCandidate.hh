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

#pragma once

#include "manta/GenomeInterval.hh"


struct SVBreakend
{
    SVBreakend() :
        isKnown(true),
        isRightOpen(true)
    {}

    // this can be set to false for one breakend when an open-ended breakend is found with unclear
    // evidence on the other side
    bool isKnown;

    // the interval is the X% confidence interval of the SV breakend, the interface allows for
    // various probabiity distributions to back this interval, but these must be accessed via
    // SVCandidate:
    GenomeInterval interval;
    bool isRightOpen;
};


struct SVCandidate
{
#if 0
    double
    breakpointProb(pos_t x, pos_t y) const;
#endif

    SVBreakend bp1;
    SVBreakend bp2;
};
