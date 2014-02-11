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

///
/// \author Chris Saunders
///

#include "manta/SVCandidateUtil.hh"



bool
isSVBelowMinSize(
    const SVCandidate& sv,
    const unsigned minSize)
{
    if (sv.bp1.interval.tid != sv.bp2.interval.tid) return false;

    const unsigned bpSize(std::abs(sv.bp1.interval.range.begin_pos() - sv.bp2.interval.range.begin_pos()));
    const unsigned insertSize(sv.insertSeq.size());

    return (std::max(bpSize,insertSize) < minSize);
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

    if (sv.bp1.interval.tid != sv.bp2.interval.tid)
    {
        return INTERTRANSLOC;
    }
    else if (SVBreakendState::isSameOrientation(sv.bp1.state,sv.bp2.state))
    {
        return INVERSION;
    }
    else if (isBp1First || isBp2First)
    {
        if (isInnies(isBp1First,sv.bp1.state,sv.bp2.state))
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
