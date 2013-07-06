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

#include "EdgeRetriever.hh"

#include "boost/foreach.hpp"

#include <cassert>

#include <iostream>

//#define DEBUG_EDGER

#ifdef DEBUG_EDGER
#include "blt_util/log.hh"
#endif



static
unsigned long
getBoundaryCount(
    const double binCount,
    const double binIndex,
    const double totalCount)
{
    return static_cast<unsigned>(std::floor((totalCount*binIndex)/binCount));
}



EdgeRetriever::
EdgeRetriever(
    const SVLocusSet& set,
    const unsigned binCount,
    const unsigned binIndex) :
    _set(set),
    _headCount(0)
{
    assert(binCount > 0);
    assert(binIndex < binCount);

    const unsigned long totalCount(_set.totalCount());
    _beginCount=(getBoundaryCount(binCount,binIndex,totalCount));
    _endCount=(getBoundaryCount(binCount,binIndex+1,totalCount));

#ifdef DEBUG_EDGER
    log_os << "EDGER: bi,bc,begin,end: "
           << binIndex << " "
           << binCount << " "
           << _beginCount << " "
           << _endCount << "\n";
#endif
}



void
EdgeRetriever::
jumpToFirstEdge()
{
    // first catch headCount up to the begin edge if required:
    while (true)
    {
        const SVLocus& locus(_set.getLocus(_edge.locusIndex));
        const unsigned locusCount(locus.totalCount());

        if ((_headCount+locusCount) > _beginCount)
        {
            while (true)
            {
                const SVLocusNode& node(locus.getNode(_edge.nodeIndex1));

                typedef SVLocusNode::edges_type::const_iterator edgeiter_t;
                edgeiter_t edgeIter(node.edges.upper_bound(_edge.nodeIndex1));
                const edgeiter_t edgeiterEnd(node.edges.end());

                for (; edgeIter != edgeiterEnd; ++edgeIter)
                {
                    const unsigned edgeCount(edgeIter->second.count + locus.getEdge(edgeIter->first,_edge.nodeIndex1).count);
                    _headCount += edgeCount;
                    if (_headCount >= _beginCount)
                    {
                        _edge.nodeIndex2 = edgeIter->first;
                        return;
                    }
                }

                _edge.nodeIndex1++;
            }
            assert(_headCount >= _beginCount);
        }
        _headCount += locusCount;
        _edge.locusIndex++;
    }
}



void
EdgeRetriever::
advanceEdge()
{
    typedef SVLocusNode::edges_type::const_iterator edgeiter_t;

    while (true)
    {
        const SVLocus& locus(_set.getLocus(_edge.locusIndex));
        while (_edge.nodeIndex1<locus.size())
        {
            const SVLocusNode& node(locus.getNode(_edge.nodeIndex1));
            edgeiter_t edgeIter(node.edges.upper_bound(_edge.nodeIndex2));
            const edgeiter_t edgeIterEnd(node.edges.end());

            for (; edgeIter != edgeIterEnd; ++edgeIter)
            {
                const unsigned edgeCount(edgeIter->second.count + locus.getEdge(edgeIter->first,_edge.nodeIndex1).count);
                _headCount += edgeCount;
                _edge.nodeIndex2 = edgeIter->first;
                return;
            }
            _edge.nodeIndex1++;
            _edge.nodeIndex2=_edge.nodeIndex1;
        }
        _edge.locusIndex++;
        _edge.nodeIndex1=0;
        _edge.nodeIndex2=0;
    }
}



bool
EdgeRetriever::
next()
{
#ifdef DEBUG_EDGER
    log_os << "EDGER: start next hc: " << _headCount << "\n";
#endif

    if (_headCount >= _endCount) return false;

    // first catch headCount up to the begin edge if required:
    if (_headCount < _beginCount)
    {
        jumpToFirstEdge();
#ifdef DEBUG_EDGER
        log_os << "EDGER: jumped hc: " << _headCount << " " << _edge  << "\n";
#endif
    }
    else
    {
        advanceEdge();
#ifdef DEBUG_EDGER
        log_os << "EDGER: advanced hc: " << _headCount << " " << _edge  << "\n";
#endif
    }

    assert(_headCount >= _beginCount);
    return (_headCount <= _endCount);
}
