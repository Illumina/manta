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
#include "EdgeOptions.hh"


/// provide an iterator over edges in a set of SV locus graphs
///
/// designed to allow parallelization of the graph processing by
/// dividing iteration into a set of bins with similar total edge
/// observation counts
///
struct EdgeRetrieverLocus : public EdgeRetriever
{
    /// \param locusIndex iterate over all edges of a specific locus
    EdgeRetrieverLocus(
        const SVLocusSet& set,
        const unsigned graphNodeMaxEdgeCount,
        const LocusEdgeOptions& opt);

    bool
    next();

private:
    void
    advanceEdge();

    LocusEdgeOptions _opt;
    bool _isInit;
};
