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
    const unsigned locusIndex) :
    EdgeRetriever(set),
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

    if(_isInit)
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
        const SVLocusNode& node(locus.getNode(_edge.nodeIndex1));
        edgeiter_t edgeIter(node.edges.lower_bound(_edge.nodeIndex2));
        const edgeiter_t edgeIterEnd(node.edges.end());

        for (; edgeIter != edgeIterEnd; ++edgeIter)
        {
            _edge.nodeIndex2 = edgeIter->first;
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
