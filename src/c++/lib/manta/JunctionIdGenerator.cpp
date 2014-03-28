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

#include "manta/JunctionIdGenerator.hh"

#include "manta/SVCandidateUtil.hh"



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


static
bool
isIntrachromBnd(
    const SVCandidate& sv)
{
    static const int intrachromTranslocThreshold(1000000);
    return (sv.centerSize() >= intrachromTranslocThreshold);
}



index_t
getExtendedSVType(
    const SVCandidate& sv)
{
    const SV_TYPE::index_t svType(getSVType(sv));

    if (svType == SV_TYPE::INTERTRANSLOC) return INTERTRANSLOC;

    if (isIntrachromBnd(sv)) return INTRATRANSLOC;

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

}


void
JunctionIdGenerator::
getId(
    const EdgeInfo& edge,
    const SVCandidate& sv,
    SVId& svId)
{
    using namespace EXTENDED_SV_TYPE;

    svId.svType=(getExtendedSVType(sv));

    svId.localId = str(_SVIdFormatter % label(svId.svType) % edge.locusIndex % edge.nodeIndex1 % edge.nodeIndex2
                       % sv.candidateIndex %  sv.assemblyAlignIndex % sv.assemblySegmentIndex );

    if      (isSVTransloc(svId.svType))
    {
        svId.mateId = svId.localId + ":1";
        svId.localId = svId.localId + ":0";
    }
    else
    {
        svId.mateId.clear();
    }
}
