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

#pragma once

#include "svgraph/SVLocus.hh"

#include <iosfwd>


struct EdgeInfo
{
    /// minimal ascii representation:
    void
    write(std::ostream& os) const;

    bool
    isSelfEdge() const
    {
        return (nodeIndex1 == nodeIndex2);
    }

    LocusIndexType locusIndex = 0;
    NodeIndexType nodeIndex1 = 0;
    NodeIndexType nodeIndex2 = 0;
};

std::ostream&
operator<<(std::ostream& os, const EdgeInfo& ei);
