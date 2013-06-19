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


#pragma once

#include "manta/SVLocus.hh"

#include "boost/foreach.hpp"


void
SVLocusNode::
mergeNode(SVLocusNode& inputNode) {
    assert(interval.tid == inputNode.interval.tid);
    assert(interval.range.is_range_intersect(inputNode.interval.range));

    interval.range.merge_range(inputNode.interval.range);
    count += inputNode.count;

    BOOST_FOREACH(edges_type::value_type& inputNodeEdgeIter, inputNode._edges)
    {
        // update local edge:
        {
            edges_type::iterator thisEdgeIter(_edges.find(inputNodeEdgeIter.first));
            if(thisEdgeIter == _edges.end())
            {
                // no self-edges:
                if(inputNodeEdgeIter.first != this)
                {
                    _edges.insert(inputNodeEdgeIter);
                }
            }
            else
            {
                thisEdgeIter->second.mergeEdge(inputNodeEdgeIter.second);
            }
        }

        // update remote inputNodeEdgeIter
        {
            edges_type& remoteEdges(inputNodeEdgeIter.first->_edges);
            edges_type::iterator oldRemoteIter(remoteEdges.find(&inputNode));
            assert(oldRemoteIter != remoteEdges.end());

            edges_type::iterator newRemoteIter(remoteEdges.find(this));

            if(newRemoteIter == remoteEdges.end())
            {
                remoteEdges.insert(std::make_pair(this,oldRemoteIter->second));
            }
            else
            {
                newRemoteIter->second.mergeEdge(oldRemoteIter->second);
            }

            remoteEdges.erase(oldRemoteIter);
        }
    }
}



void
SVLocusNode::
clearEdges()
{
    BOOST_FOREACH(edges_type::value_type& edgeIter, _edges)
    {
        edges_type& remoteEdges(edgeIter.first->_edges);
        edges_type::iterator thisRemoteIter(remoteEdges.find(this));
        assert(thisRemoteIter != remoteEdges.end());

        remoteEdges.erase(thisRemoteIter);
    }
    _edges.clear();
}
