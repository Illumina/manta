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

#pragma once

#include "svgraph/EdgeInfo.hh"
#include "svgraph/SVLocusSet.hh"


/// provide an iterator over edges in a set of SV locus graphs
///
struct EdgeRetriever
{
    EdgeRetriever(
        const SVLocusSet& set,
        const unsigned graphNodeMaxEdgeCount) :
        _set(set),
        _graphNodeMaxEdgeCount(graphNodeMaxEdgeCount)
    {}

    virtual
    ~EdgeRetriever()
    {}

    virtual
    bool
    next() = 0;

    const EdgeInfo&
    getEdge() const
    {
        return _edge;
    }

protected:
    const SVLocusSet& _set;
    const unsigned _graphNodeMaxEdgeCount;
    EdgeInfo _edge;
};
