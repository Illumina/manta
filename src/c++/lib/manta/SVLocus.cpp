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


#include "common/Exceptions.hh"
#include "manta/SVLocus.hh"

#include "boost/foreach.hpp"

#include <iostream>



std::ostream&
operator<<(std::ostream& os, const GenomeInterval& gi)
{
    os << "GenomeIntreval: " << gi.tid << ":" << gi.range;
    return os;
}


std::ostream&
operator<<(std::ostream& os, const SVLocusEdge& edge)
{
    os << "Edgecount: " << edge.count;
    return os;
}



void
SVLocusNode::
mergeNode(SVLocusNode& inputNode) {

    using namespace illumina::common;

    if ((interval.tid != inputNode.interval.tid) ||
        (! interval.range.is_range_intersect(inputNode.interval.range))) {
        std::ostringstream oss;
        oss << "ERROR: Attempting to merge non-intersecting nodes\n"
            << "\tNode1: " << *this
            << "\tNode2: " << inputNode;
        BOOST_THROW_EXCEPTION(PreConditionException(oss.str()));
    }

    interval.range.merge_range(inputNode.interval.range);
    count += inputNode.count;

    BOOST_FOREACH(edges_type::value_type& inputNodeEdgeIter, inputNode._edges)
    {
        // update local edge:
        {
            edges_type::iterator thisEdgeIter(_edges.find(inputNodeEdgeIter.first));
            if (thisEdgeIter == _edges.end())
            {
                // no self-edges:
                if (inputNodeEdgeIter.first != this)
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

            if (newRemoteIter == remoteEdges.end())
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
    using namespace illumina::common;

    BOOST_FOREACH(edges_type::value_type& edgeIter, _edges)
    {
        edges_type& remoteEdges(edgeIter.first->_edges);
        edges_type::iterator thisRemoteIter(remoteEdges.find(this));
        if(thisRemoteIter == remoteEdges.end())
        {
            std::ostringstream oss;
            oss << "ERROR: no return edge on remote node.\n";
            oss << "\tlocal_node: " << *this;
            oss << "\tremote_ndoe: " << *(edgeIter.first);
            BOOST_THROW_EXCEPTION(PreConditionException(oss.str()));
        }

        remoteEdges.erase(thisRemoteIter);
    }
    _edges.clear();
}



std::ostream&
operator<<(std::ostream& os, const SVLocusNode& node)
{
    os << "LocusNode: " << &node << " count: " << node.count << " " << node.interval << "\n";

    BOOST_FOREACH(const SVLocusNode::edges_type::value_type& edgeIter, node._edges)
    {
        os << "\tEdgeTo: " << edgeIter.first << " " << edgeIter.second << "\n";
    }
    return os;
}



void
SVLocus::
checkState() const
{
    assert(_graph.size() == _smap.size());

    BOOST_FOREACH(const SVLocusNode* nodePtr, *this)
    {
        assert(1 == _smap.count(nodePtr));
    }
}


std::ostream&
operator<<(std::ostream& os, const SVLocus& locus)
{
    os << "LOCUS_BEGIN\n";
    BOOST_FOREACH(const SVLocusNode* nodePtr, locus)
    {
        os << *nodePtr;
    }
    os << "LOCUS_END\n";
    return os;
}
