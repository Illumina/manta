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

#pragma once

#include "manta/SVLocusSet.hh"

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



struct EdgeRetriever
{
    EdgeRetriever(
        const unsigned binIndex,
        const unsigned binCount,
        const SVLocusSet& set);

    bool
    next();

    const EdgeInfo&
    getEdge() const
    {
        return _edge;
    }

private:
    void
    jumpToFirstEdge();

    void
    advanceEdge();

    const SVLocusSet& _set;
    unsigned long _beginCount;
    unsigned long _endCount;

    unsigned long _headCount;
    EdgeInfo _edge;
};

