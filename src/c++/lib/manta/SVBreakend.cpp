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

///
/// \author Chris Saunders
///

#include "manta/SVBreakend.hh"

#include <iostream>


std::ostream&
operator<<(
    std::ostream& os,
    const SVBreakendLowResEvidence& sce)
{
    os << "SVBreakendLowResEvidence:";
    for (unsigned i(0); i<SVBreakendLowResEvidence::SIZE; ++i)
    {
        os << " " << SVBreakendLowResEvidence::label(i) << ": " << sce.getVal(i);
    }

    return os;
}


std::ostream&
operator<<(
    std::ostream& os,
    const SVBreakend& svb)
{
    os << "Breakend:"
       << " " << svb.interval
       << " " << SVBreakendState::label(svb.state) << "\n"
       << "\t" << svb.lowresEvidence << "\n";
    return os;
}
