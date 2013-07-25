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

#include "EdgeRetriever.hh"


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
        const unsigned locusIndex);

    bool
    next();

private:
    void
    advanceEdge();

    unsigned _locusIndex;
    bool _isInit;
};
