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

#include "manta/SVCandidate.hh"

#include <iostream>



std::ostream&
operator<<(std::ostream& os, const SVBreakend& svb)
{
    os << "Breakend:"
       << " " << svb.interval
       << " " << SVBreakendState::label(svb.state)
       << " readCount: " << svb.readCount
       << " pairCount: " << svb.pairCount;
    return os;
}



std::ostream&
operator<<(std::ostream& os, const SVCandidate& svc)
{
    os << "SVCandidate:\n"
       << "\t" << svc.bp1 << "\n"
       << "\t" << svc.bp2 << "\n";
    return os;
}
