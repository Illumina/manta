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

#include "EdgeRetrieverLocus.hh"

#include "boost/foreach.hpp"

#include <cassert>

#include <iostream>

//#define DEBUG_EDGER

#ifdef DEBUG_EDGER
#include "blt_util/log.hh"
#endif



EdgeRetrieverLocus::
EdgeRetrieverLocus(
    const SVLocusSet& set,
    const unsigned graphNodeMaxEdgeCount,
    const unsigned locusIndex) :
    EdgeRetriever(set, graphNodeMaxEdgeCount),
    _locusIndex(locusIndex),
    _isInit(false)
{
    assert(_locusIndex<set.size());
}



void
EdgeRetrieverLocus::
advanceEdge()
{
    typedef SVLocusNode::edges_type::const_iterator edgeiter_t;

    if (_isInit)
    {
        _edge.nodeIndex2++;
    }
    else
    {
        _edge.locusIndex=_locusIndex;
        _edge.nodeIndex1=0;
        _edge.nodeIndex2=0;
        _isInit=true;
    }

    const SVLocus& locus(_set.getLocus(_edge.locusIndex));
    while (_edge.nodeIndex1<locus.size())
    {
        const SVLocusNode& node1(locus.getNode(_edge.nodeIndex1));
        const bool isEdgeFilterNode1((_graphNodeMaxEdgeCount>0) && node1.edges.size()>_graphNodeMaxEdgeCount);
        edgeiter_t edgeIter(node1.edges.lower_bound(_edge.nodeIndex2));
        const edgeiter_t edgeIterEnd(node1.edges.end());

        for (; edgeIter != edgeIterEnd; ++edgeIter)
        {
            _edge.nodeIndex2 = edgeIter->first;

            // check whether this is a noise edge that we skip:
            if (isEdgeFilterNode1)
            {
                const SVLocusNode& node2(locus.getNode(_edge.nodeIndex2));
                const bool isEdgeFilterNode2(node2.edges.size()>_graphNodeMaxEdgeCount);
                if (isEdgeFilterNode2) continue;
            }
            return;
        }
        _edge.nodeIndex1++;
        _edge.nodeIndex2=_edge.nodeIndex1;
    }

    _edge.locusIndex++;
}



bool
EdgeRetrieverLocus::
next()
{
#ifdef DEBUG_EDGER
    log_os << "EDGERL: start\n";
#endif

    if (_edge.locusIndex > _locusIndex) return false;
    advanceEdge();

    return (_edge.locusIndex == _locusIndex);
}
