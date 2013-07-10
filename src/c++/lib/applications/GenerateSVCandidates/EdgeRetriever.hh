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

#include "svgraph/EdgeInfo.hh"
#include "svgraph/SVLocusSet.hh"


/// provide an iterator over edges in a set of SV locus graphs
///
/// designed to allow parallelization of the graph processing by
/// dividing iteration into a set of bins with similar total edge
/// observation counts
///
struct EdgeRetriever
{
    /// \param binCount total number of parallel bins, must be 1 or greater
    /// \param binIndex parallel bin id, must be less than binCount
    EdgeRetriever(
        const SVLocusSet& set,
        const unsigned binCount = 1,
        const unsigned binIndex = 0);

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
