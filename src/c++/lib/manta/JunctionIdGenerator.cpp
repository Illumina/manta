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
