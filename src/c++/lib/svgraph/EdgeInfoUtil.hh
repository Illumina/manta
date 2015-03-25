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


inline
bool
isBidirectionalEdge(
    const SVLocusSet& cset,
    const EdgeInfo& edge)
{
    const unsigned minEdgeCount(cset.getMinMergeEdgeCount());

    const SVLocus& locus(cset.getLocus(edge.locusIndex));

    return ((locus.getEdge(edge.nodeIndex1,edge.nodeIndex2).getCount() >= minEdgeCount) &&
            (locus.getEdge(edge.nodeIndex2,edge.nodeIndex1).getCount() >= minEdgeCount));
}


/// determine if this is a self-edge of a node with no other (bidirectional-pass) edges:
bool
testIsolatedEdge(
    const SVLocusSet& cset,
    const EdgeInfo& edge);
