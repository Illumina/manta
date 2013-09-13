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

#include <cassert>

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
       << "\tisImprecise?: " << svc.isImprecise() << "\n"
       << "\tcandidate:assembly index: " << svc.candidateIndex << ":" << svc.assemblyIndex << "\n";
    if (! svc.isImprecise())
    {
       os << "\tAlignment: " << svc.insertAlignment << "\n"
          << "\tBreakendInsertSeq: " << svc.insertSeq << "\n";
    }
    os << "\t" << svc.bp1 << "\n"
       << "\t" << svc.bp2 << "\n";
    return os;
}



SV_TYPE::index_t
getSVType(const SVCandidate& sv)
{
    using namespace SV_TYPE;

    // remove failed local assemblies first:
    if ((sv.bp1.state == SVBreakendState::UNKNOWN) || (sv.bp2.state == SVBreakendState::UNKNOWN))
    {
        return UNKNOWN;
    }

    const bool isBp1First(sv.bp1.interval.range.begin_pos() <= sv.bp2.interval.range.begin_pos());
    const bool isBp2First(sv.bp2.interval.range.begin_pos() <= sv.bp1.interval.range.begin_pos());

    // true for insertions:
    const bool isBpEqual(sv.bp1.interval.range == sv.bp2.interval.range);

    assert(! (isBp1First && isBp2First));

    if (sv.bp1.interval.tid != sv.bp2.interval.tid)
    {
        return INTERTRANSLOC;
    }
    else if (SVBreakendState::isSameOrientation(sv.bp1.state,sv.bp2.state))
    {
        return INVERSION;
    }
    else if(isBpEqual || isBp1First || isBp2First)
    {
        if (isBpEqual || isInnies(isBp1First,sv.bp1.state,sv.bp2.state))
        {
            return INDEL;
        }
        else if (isOutties(isBp1First,sv.bp1.state,sv.bp2.state))
        {
            return TANDUP;
        }
    }

    return UNKNOWN;
}
