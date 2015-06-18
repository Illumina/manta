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
/// \author Chris Saunders
///

#include "manta/SVCandidateUtil.hh"



bool
isSVBelowMinSize(
    const SVCandidate& sv,
    const unsigned minSize)
{
    if (sv.bp1.interval.tid != sv.bp2.interval.tid) return false;

    const unsigned bpSize(std::abs(sv.bp1.interval.range.center_pos() - sv.bp2.interval.range.center_pos())-1);
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



namespace EXTENDED_SV_TYPE
{

/// is an indel classified as insert or delete?
static
index_t
classifyIndel(
    const SVCandidate& sv)
{
    const bool isBp1First(sv.bp1.interval.range.begin_pos()<=sv.bp2.interval.range.begin_pos());

    const SVBreakend& bpA(isBp1First ? sv.bp1 : sv.bp2);
    const SVBreakend& bpB(isBp1First ? sv.bp2 : sv.bp1);

    const unsigned deleteSize(bpB.interval.range.begin_pos() - bpA.interval.range.begin_pos());
    const unsigned insertSize(sv.insertSeq.size());

    return ((deleteSize >= insertSize) ? DELETE : INSERT);
}

}



EXTENDED_SV_TYPE::index_t
getExtendedSVType(
    const SVCandidate& sv,
    const bool isForceIntraChromBnd)
{
    using namespace EXTENDED_SV_TYPE;

    const SV_TYPE::index_t svType(getSVType(sv));

    if (svType == SV_TYPE::INTERTRANSLOC) return INTERTRANSLOC;

    if (isForceIntraChromBnd) return INTRATRANSLOC;

    switch (svType)
    {
    case SV_TYPE::INVERSION:
        return INVERSION;
    case SV_TYPE::TANDUP:
        return TANDUP;
    case SV_TYPE::INDEL:
        return classifyIndel(sv);
    default:
        return UNKNOWN;
    }
}
