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

#include "svgraph/EdgeInfoUtil.hh"



bool
testIsolatedEdge(
    const SVLocusSet& cset,
    const EdgeInfo& edge)
{
    if (edge.nodeIndex1 != edge.nodeIndex2) return false;

    const SVLocus& locus(cset.getLocus(edge.locusIndex));
#if 0
    // simple criteria -- make sure there are no other nodes in locus
    return (locus.size() == 1);
#endif

    // search to check to see if any bidirectional edges extend from this node (other than the self-edge):
    const SVLocusNode& node1(locus.getNode(edge.nodeIndex1));
    const SVLocusEdgeManager node1Manager(node1.getEdgeManager());

    typedef SVLocusEdgesType::const_iterator edgeiter_t;
    edgeiter_t edgeIter(node1Manager.getMap().begin());
    const edgeiter_t edgeIterEnd(node1Manager.getMap().end());

    EdgeInfo testEdge(edge);

    for (; edgeIter != edgeIterEnd; ++edgeIter)
    {
        testEdge.nodeIndex2 = (edgeIter->first);
        if (testEdge.nodeIndex1 == testEdge.nodeIndex2) continue;
        if (isBidirectionalEdge(cset, testEdge)) return false;
    }

    return true;
}
