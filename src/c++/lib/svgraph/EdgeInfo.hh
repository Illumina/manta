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

#pragma once

#include "svgraph/SVLocus.hh"

#include <iosfwd>


struct EdgeInfo
{
    EdgeInfo() :
        locusIndex(0),
        nodeIndex1(0),
        nodeIndex2(0)
    {}

    LocusIndexType locusIndex;
    NodeIndexType nodeIndex1;
    NodeIndexType nodeIndex2;
};

std::ostream&
operator<<(std::ostream& os, const EdgeInfo& ei);
