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

#include "common/Exceptions.hh"
#include "manta/SVLocus.hh"

#include <iostream>

#ifdef DEBUG_SVL
#include "blt_util/log.hh"
#endif



std::ostream&
operator<<(std::ostream& os, const GenomeInterval& gi)
{
    os << "GenomeInterval: " << gi.tid << ":" << gi.range;
    return os;
}


std::ostream&
operator<<(std::ostream& os, const SVLocusEdge& edge)
{
    os << "Edgecount: " << edge.count;
    return os;
}


void
SVLocus::
mergeNode(
        const NodeIndexType fromIndex,
        const NodeIndexType toIndex)
{
    using namespace illumina::common;

#ifdef DEBUG_SVL
    log_os << "mergeNode from: " << fromIndex << " to: " << toIndex << " size: " << size() << "\n";
#endif

    SVLocusNode& fromNode(getNode(fromIndex));
    SVLocusNode& toNode(getNode(toIndex));

    if(fromNode.interval.tid != toNode.interval.tid)
    {
        std::ostringstream oss;
        oss << "ERROR: Attempting to merge nodes on different chromosomes\n"
            << "\tNode1: " << fromNode
            << "\tNode2: " << toNode;
        BOOST_THROW_EXCEPTION(PreConditionException(oss.str()));
    }

    notifyDelete(toIndex);

    toNode.interval.range.merge_range(fromNode.interval.range);
    toNode.count += fromNode.count;

    notifyAdd(toIndex);

    BOOST_FOREACH(edges_type::value_type& fromNodeEdgeIter, fromNode)
    {

#ifdef DEBUG_SVL
        log_os << "mergeNode: handle fromEdge: " << _index << ":" << fromNodeEdgeIter.first << "\n";
#endif

        // update local edge:
        {
            // this node and input node could both have an edge to the same node already
            // -- check that case here:
            edges_type::iterator toNodeEdgeIter(toNode.edges.find(fromNodeEdgeIter.first));
            if (toNodeEdgeIter == toNode.edges.end())
            {
#ifdef DEBUG_SVL
                log_os << "mergeNode: fromEdge is not in toNode\n";
#endif
                // no self-edges:
                if (fromNodeEdgeIter.first != toIndex)
                {
#ifdef DEBUG_SVL
                    log_os << "mergeNode: toNode add " << _index << ":" << fromNodeEdgeIter.first << "\n";
#endif

                    toNode.edges.insert(fromNodeEdgeIter);
                }
            }
            else
            {
#ifdef DEBUG_SVL
                log_os << "mergeNode: fromEdge is already in toNode\n";
#endif
                // this node does contain a link to the remote node already:
                toNodeEdgeIter->second.mergeEdge(fromNodeEdgeIter.second);
            }
        }

        // update remote inputNodeEdgeIter
        {
            SVLocusNode& remoteNode(getNode(fromNodeEdgeIter.first));
            edges_type& remoteEdges(remoteNode.edges);
            edges_type::iterator oldRemoteIter(remoteEdges.find(fromIndex));
            if(oldRemoteIter == remoteEdges.end())
            {
                std::ostringstream oss;
                oss << "ERROR: Can't find return edge to node index: " << _index << ":" << fromIndex << " in remote node index: " << _index << ":" << fromNodeEdgeIter.first << "\n"
                    << "\tlocal_node: " << fromNode
                    << "\tremote_node: " << remoteNode;
                BOOST_THROW_EXCEPTION(PreConditionException(oss.str()));
            }

            // the remote node could contain a link to toIndex already, check that here:
            edges_type::iterator newRemoteIter(remoteEdges.find(toIndex));
            if (newRemoteIter == remoteEdges.end())
            {
#ifdef DEBUG_SVL
                log_os << "mergeNode: fromRemote does not point to toIndex\n";
#endif
                remoteEdges.insert(std::make_pair(toIndex,oldRemoteIter->second));
            }
            else
            {
#ifdef DEBUG_SVL
                log_os << "mergeNode: fromRemote already points to toIndex\n";
#endif
                newRemoteIter->second.mergeEdge(oldRemoteIter->second);
            }
        }
    }

    clearNodeEdges(fromIndex);
}



void
SVLocus::
clearNodeEdges(NodeIndexType nodePtr)
{
    using namespace illumina::common;

    SVLocusNode& node(getNode(nodePtr));
    BOOST_FOREACH(edges_type::value_type& edgeIter, node)
    {
        SVLocusNode& remoteNode(getNode(edgeIter.first));
        edges_type& remoteEdges(remoteNode.edges);
        edges_type::iterator thisRemoteIter(remoteEdges.find(nodePtr));
        if(thisRemoteIter == remoteEdges.end())
        {
            std::ostringstream oss;
            oss << "ERROR: no return edge on remote node.\n";
            oss << "\tlocal_node: " << node;
            oss << "\tremote_node: " << remoteNode;
            BOOST_THROW_EXCEPTION(PreConditionException(oss.str()));
        }

        remoteEdges.erase(thisRemoteIter);
    }
    node.edges.clear();
}



void
SVLocus::
eraseNode(const NodeIndexType nodePtr)
{
    using namespace illumina::common;

    if (nodePtr >= _graph.size()) return;

    clearNodeEdges(nodePtr);

    NodeIndexType fromPtr(_graph.size()-1);

#ifdef DEBUG_SVL
    log_os << "eraseNode: from: " << _index << ":" << fromPtr << " to: " << _index << ":" << nodePtr << "\n";
#endif

    if(fromPtr != nodePtr)
    {
        // reassign fromNode's remote edges before shifting its address:
        //
        SVLocusNode& fromNode(getNode(fromPtr));
        BOOST_FOREACH(edges_type::value_type& edgeIter, fromNode)
        {
            SVLocusNode& remoteNode(getNode(edgeIter.first));
            edges_type& remoteEdges(remoteNode.edges);
            edges_type::iterator thisRemoteIter(remoteEdges.find(fromPtr));
            if(thisRemoteIter == remoteEdges.end())
            {
                std::ostringstream oss;
                oss << "ERROR: eraseNode: no return edge on remote node.\n";
                oss << "\tlocal_node: " << fromNode;
                oss << "\tremote_node: " << remoteNode;
                BOOST_THROW_EXCEPTION(PreConditionException(oss.str()));
            }

            remoteEdges.insert(std::make_pair(nodePtr,thisRemoteIter->second));
            remoteEdges.erase(thisRemoteIter);
        }

        notifyDelete(nodePtr);
        _graph[nodePtr] = _graph[fromPtr];
        notifyAdd(nodePtr);
    }
    notifyDelete(fromPtr);
    _graph.resize(fromPtr);
}



std::ostream&
operator<<(std::ostream& os, const SVLocusNode& node)
{
    os << "LocusNode: count: " << node.count << " " << node.interval << "\n";

    BOOST_FOREACH(const SVLocusNode::edges_type::value_type& edgeIter, node)
    {
        os << "\tEdgeTo: " << edgeIter.first << " " << edgeIter.second << "\n";
    }
    return os;
}



void
SVLocus::
checkState() const
{
    using namespace illumina::common;

    const unsigned nodeSize(size());
    for(unsigned nodeIndex(0);nodeIndex<nodeSize;++nodeIndex)
    {
        const SVLocusNode& node(getNode(nodeIndex));

        BOOST_FOREACH(const edges_type::value_type& edgeIter, node)
        {
            // check that that every edge has a return path:
            const SVLocusNode& edgeNode(getNode(edgeIter.first));
            SVLocusNode::const_iterator iter(edgeNode.edges.find(nodeIndex));
            if(edgeNode.edges.end() == iter)
            {
                std::ostringstream oss;
                oss << "ERROR: no return edge on remote node.\n";
                oss << "\tlocal_node: " << _index << ":" << nodeIndex << " " << node;
                oss << "\tremote_node: " << _index << ":" << edgeIter.first << " " << edgeNode;
                BOOST_THROW_EXCEPTION(PreConditionException(oss.str()));
            }
        }
    }
}



std::ostream&
operator<<(std::ostream& os, const SVLocus& locus)
{
    os << "LOCUS BEGIN INDEX " << locus.getIndex() << "\n";
    unsigned locusIndex(0);
    BOOST_FOREACH(const SVLocusNode& node, locus)
    {
        os << "LocusIndex: " << locusIndex << " " << node;
        locusIndex++;
    }
    os << "LOCUS END INDEX " << locus.getIndex() << "\n";
    return os;
}
