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
    using namespace SVEvidenceType;

    os << "SVBreakendLowResEvidence:";
    for (int i(0); i<SIZE; ++i)
    {
        os << " " << label(i) << ": " << sce.getVal(i);
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
