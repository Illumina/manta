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

#include "EdgeRetriever.hh"


/// provide an iterator over edges in a set of SV locus graphs
///
/// designed to allow parallelization of the graph processing by
/// dividing iteration into a set of bins with similar total edge
/// observation counts
///
/// the contents of the bins are designed to be distributed evenly over the sequence of edges
///
struct EdgeRetrieverJumpBin : public EdgeRetriever
{
    /// \param[in] graphNodeMaxEdgeCount filtration parameter for skipping edges from highly connected nodes (set to zero to disable)
    /// \param[in] binCount total number of parallel bins, must be 1 or greater
    /// \param[in] binIndex parallel bin id, must be less than binCount
    EdgeRetrieverJumpBin(
        const SVLocusSet& set,
        const unsigned graphNodeMaxEdgeCount,
        const unsigned binCount,
        const unsigned binIndex);

    bool
    next();

private:
    void
    advanceEdge();

    typedef unsigned long count_t;

    unsigned _binCount;
    unsigned _binIndex;
    count_t _edgeIndex;

    // additional 'load balancing' structures:
    count_t _avgBinTotalCount;

   std::vector<count_t> _binTotalCount;
};
