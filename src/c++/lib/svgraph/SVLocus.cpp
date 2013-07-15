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
#include "svgraph/SVLocus.hh"

#include <iostream>
#include <stack>

#ifdef DEBUG_SVL
#include "blt_util/log.hh"
#endif



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

    if (fromNode.interval.tid != toNode.interval.tid)
    {
        std::ostringstream oss;
        oss << "ERROR: Attempting to merge nodes on different chromosomes\n"
            << "\tNode1: " << fromNode
            << "\tNode2: " << toNode;
        BOOST_THROW_EXCEPTION(LogicException(oss.str()));
    }

    notifyDelete(toIndex);

    toNode.interval.range.merge_range(fromNode.interval.range);
    if     ((toNode.count==0) && (fromNode.count!=0))
    {
        toNode.evidenceRange = fromNode.evidenceRange;
    }
    else if ((fromNode.count==0) && (toNode.count!=0))
    {
        // pass (keep toNode value as is
    }
    else
    {
        toNode.evidenceRange.merge_range(fromNode.evidenceRange);
    }
    toNode.count += fromNode.count;

    notifyAdd(toIndex);

    BOOST_FOREACH(edges_type::value_type& fromNodeEdgeIter, fromNode)
    {
#ifdef DEBUG_SVL
        // is this edge between the to and from nodes?
        const bool isToFromEdge(fromNodeEdgeIter.first == toIndex);
       log_os << "mergeNode: handle fromEdge: " << _index << ":" << fromNodeEdgeIter.first << " isToFromEdge: " << isToFromEdge << "\n";
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
                log_os << "mergeNode: toNode add " << _index << ":" << fromNodeEdgeIter.first << "\n";
#endif

                toNode.edges.insert(fromNodeEdgeIter);
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
            if (oldRemoteIter == remoteEdges.end())
            {
                std::ostringstream oss;
                oss << "ERROR: Can't find return edge to node index: " << _index << ":" << fromIndex << " in remote node index: " << _index << ":" << fromNodeEdgeIter.first << "\n"
                    << "\tlocal_node: " << fromNode
                    << "\tremote_node: " << remoteNode;
                BOOST_THROW_EXCEPTION(LogicException(oss.str()));
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
getEdgeException(
    const NodeIndexType fromIndex,
    const NodeIndexType toIndex) const
{
    using namespace illumina::common;

    std::ostringstream oss;
    oss << "ERROR: SVLocus::getEdge() no edge exists\n";
    oss << "\tfrom_node: " << _index << ":" << fromIndex << " " << getNode(fromIndex);
    oss << "\tto_node: " << _index << ":" << toIndex << " " << getNode(toIndex);
    BOOST_THROW_EXCEPTION(LogicException(oss.str()));
}



bool
SVLocus::
isNoiseNode(
    const unsigned minMergeEdgeCount,
    const NodeIndexType nodeIndex) const
{
    const SVLocusNode& node(getNode(nodeIndex));
    BOOST_FOREACH(const SVLocusNode::edges_type::value_type& edge, node)
    {
        if (edge.second.count >= minMergeEdgeCount) return false;
        if (getEdge(edge.first,nodeIndex).count >= minMergeEdgeCount) return false;
    }
    return true;
}



unsigned
SVLocus::
cleanNodeCore(
    const unsigned minMergeEdgeCount,
    const NodeIndexType nodeIndex,
    std::set<NodeIndexType>& emptyNodes)
{
#ifdef DEBUG_SVL
    log_os << "cleanNodeCore nodeAddy: " << _index << ":" << nodeIndex << "\n";
#endif

    unsigned totalCleaned(0);
    SVLocusNode& queryNode(getNode(nodeIndex));

    std::vector<NodeIndexType> eraseEdges;
    BOOST_FOREACH(edges_type::value_type& edgeIter, queryNode)
    {
        if (0 != edgeIter.second.count)
        {
            if (edgeIter.second.count < minMergeEdgeCount)
            {
                // clean criteria met -- go ahead and erase edge count:
                assert(queryNode.count>=edgeIter.second.count);
                totalCleaned += edgeIter.second.count;
                queryNode.count -= edgeIter.second.count;
                edgeIter.second.count = 0;
            }
        }

        if (0 == edgeIter.second.count)
        {
            // if the out edge count is zero, see if the in-edge count is also zero --
            // if so, erase edge
            //
            const SVLocusEdge& fromRemoteEdge(getEdge(edgeIter.first,nodeIndex));
            if (0 == fromRemoteEdge.count)
            {
                eraseEdges.push_back(edgeIter.first);

                // also check to see if the remote node will be empty after
                // this edge deletion:
                const SVLocusNode& remoteNode(getNode(edgeIter.first));
                if ((0 == remoteNode.count) &&
                    (1 == remoteNode.edges.size()))
                {
                    emptyNodes.insert(edgeIter.first);
                }
            }
        }
    }

    // delete empty edges:
    BOOST_FOREACH(const NodeIndexType toIndex, eraseEdges)
    {
        clearEdgePair(nodeIndex,toIndex);
    }

    // if true add the target node to the erase list:
    if ((0 == queryNode.edges.size()) && (0 == queryNode.count))
    {
        emptyNodes.insert(nodeIndex);
    }

#ifdef DEBUG_SVL
    log_os << "cleanNodeCore emptyEdges:\n";
    BOOST_FOREACH(const NodeIndexType toIndex, eraseEdges)
    {
        log_os << "\tedge: " << _index << ":" << nodeIndex << "->" << _index << ":" << toIndex << "\n";
    }

    log_os << "cleanNodeCore emptyNodes\n";
    BOOST_FOREACH(const NodeIndexType nodeIndex2, emptyNodes)
    {
        log_os << "\tnodeAddy: " << _index << ":" << nodeIndex2 << "\n";
    }

    log_os << "totalCleaned: " << totalCleaned << "\n";
#endif

    return totalCleaned;
}



unsigned
SVLocus::
cleanNode(
    const unsigned minMergeEdgeCount,
    const NodeIndexType nodeIndex)
{
    std::set<NodeIndexType> emptyNodes;
    const unsigned totalCleaned(cleanNodeCore(minMergeEdgeCount,nodeIndex,emptyNodes));
    eraseNodes(emptyNodes);
    return totalCleaned;
}



unsigned
SVLocus::
clean(const unsigned minMergeEdgeCount)
{
    std::set<NodeIndexType> emptyNodes;
    unsigned totalCleaned(0);

    const unsigned nodeSize(size());
    for (unsigned nodeIndex(0); nodeIndex<nodeSize; ++nodeIndex)
    {
        totalCleaned += cleanNodeCore(minMergeEdgeCount,nodeIndex,emptyNodes);
    }
    eraseNodes(emptyNodes);
    return totalCleaned;
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
        if (thisRemoteIter == remoteEdges.end())
        {
            std::ostringstream oss;
            oss << "ERROR: no return edge on remote node.\n";
            oss << "\tlocal_node: " << node;
            oss << "\tremote_node: " << remoteNode;
            BOOST_THROW_EXCEPTION(LogicException(oss.str()));
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
    log_os << "eraseNode: " << _index << ":" << nodePtr << " transfer_in: " << _index << ":" << fromPtr << " \n";
#endif

    if (fromPtr != nodePtr)
    {
        // reassign fromNode's remote edges before shifting its address:
        //
        SVLocusNode& fromNode(getNode(fromPtr));
        BOOST_FOREACH(edges_type::value_type& edgeIter, fromNode)
        {
            SVLocusNode& remoteNode(getNode(edgeIter.first));
            edges_type& remoteEdges(remoteNode.edges);
            remoteEdges.insert(std::make_pair(nodePtr,getEdge(edgeIter.first,fromPtr)));
            remoteEdges.erase(fromPtr);
        }

        notifyDelete(nodePtr);
        _graph[nodePtr] = _graph[fromPtr];
        notifyAdd(nodePtr);
    }
    notifyDelete(fromPtr);
    _graph.resize(fromPtr);
}



void
SVLocus::
eraseNodes(const std::set<NodeIndexType>& nodes)
{

    if (nodes.empty()) return;

    if (size() == nodes.size())
    {
        // if the whole locus is being erased, this is more efficient:
        clear();
        return;
    }

    // partial deletion must be done in descending order:
    BOOST_REVERSE_FOREACH(const NodeIndexType nodeIndex, nodes)
    {
        eraseNode(nodeIndex);
    }
}



std::ostream&
operator<<(std::ostream& os, const SVLocusNode& node)
{
    os << "LocusNode: count: " << node.count << " " << node.interval
       << " n_edges: " << node.size()
       << " out_count: " << node.outCount()
       << " evidence: " << node.evidenceRange
       << "\n";

    BOOST_FOREACH(const SVLocusNode::edges_type::value_type& edgeIter, node)
    {
        os << "\tEdgeTo: " << edgeIter.first
           << " out_count: " << edgeIter.second.count
           << "\n";
    }
    return os;
}



unsigned
SVLocus::
getNodeInCount(
    const LocusIndexType nodeIndex) const
{
    const SVLocusNode& node(getNode(nodeIndex));

    unsigned sum(0);
    BOOST_FOREACH(const SVLocusNode::edges_type::value_type& edgeIter, node)
    {
        sum += getEdge(edgeIter.first,nodeIndex).count;
    }
    return sum;
}




void
SVLocus::
dumpNode(
    std::ostream& os,
    const LocusIndexType nodeIndex) const
{
    const SVLocusNode& node(getNode(nodeIndex));
    os << "LocusNode: count: " << node.count << " " << node.interval
       << " n_edges: " << node.size()
       << " out_count: " << node.outCount()
       << " in_count: " << getNodeInCount(nodeIndex)
       << " evidence: " << node.evidenceRange
       << "\n";

    BOOST_FOREACH(const SVLocusNode::edges_type::value_type& edgeIter, node)
    {
        os << "\tEdgeTo: " << edgeIter.first
           << " out_count: " << edgeIter.second.count
           << " in_count: " << getEdge(edgeIter.first,nodeIndex).count << "\n";
    }
}



void
SVLocus::
findConnected(
    const NodeIndexType startIndex,
    std::set<NodeIndexType>& connected) const
{
    connected.clear();

    std::stack<NodeIndexType> nodeStack;
    nodeStack.push(startIndex);

    while (! nodeStack.empty())
    {
        connected.insert(nodeStack.top());
        const SVLocusNode& node(getNode(nodeStack.top()));
        nodeStack.pop();
        BOOST_FOREACH(const edges_type::value_type& edgeIter, node)
        {
            if (! connected.count(edgeIter.first)) nodeStack.push(edgeIter.first);
        }
    }
}



void
SVLocus::
checkState(const bool isCheckConnected) const
{
    using namespace illumina::common;

    const unsigned nodeSize(size());
    for (unsigned nodeIndex(0); nodeIndex<nodeSize; ++nodeIndex)
    {
        const SVLocusNode& node(getNode(nodeIndex));

        BOOST_FOREACH(const edges_type::value_type& edgeIter, node)
        {
            // check that that every edge has a return path:
            getEdge(edgeIter.first,nodeIndex);
        }
    }

    if (0 == nodeSize) return;

    if (! isCheckConnected) return;

    // check that every locus in the graph is connected:
    std::set<NodeIndexType> connected;
    findConnected(0,connected);

    if (nodeSize != connected.size())
    {
        std::ostringstream oss;
        oss << "ERROR: SVLocus contains unconnected components, LocusIndex: " << _index << "\n";
        BOOST_THROW_EXCEPTION(LogicException(oss.str()));
    }
}



std::ostream&
operator<<(std::ostream& os, const SVLocus& locus)
{
    os << "LOCUS BEGIN INDEX " << locus.getIndex() << "\n";
    const unsigned nodeCount(locus.size());
    for (unsigned nodeIndex(0); nodeIndex<nodeCount; ++nodeIndex)
    {
        os << "NodeIndex: " << nodeIndex << " ";
        locus.dumpNode(os,nodeIndex);
    }
    os << "LOCUS END INDEX " << locus.getIndex() << "\n";
    return os;
}
