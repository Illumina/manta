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
// <https://github.com/sequencing/licenses/>
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
    os << "Edgecount: " << edge.getCount() << " isCountExact?: " << edge.isCountExact();
    return os;
}



void
SVLocusNode::
getEdgeException(
    const NodeIndexType toIndex) const
{
    using namespace illumina::common;

    std::ostringstream oss;
    oss << "ERROR: SVLocusNode::getEdge() no edge exists\n";
    oss << "\tfrom node: " << (*this) << "\n";
    oss << "\tto_node index: " << toIndex << "\n";
    BOOST_THROW_EXCEPTION(LogicException(oss.str()));
}



void
SVLocus::
nodeHurl(const NodeIndexType nodePtr) const
{
    using namespace illumina::common;

    std::ostringstream oss;
    oss << "ERROR: Attempting to access node: " << _index << ":" << nodePtr << " in locus with size: " << size() << "\n";
    BOOST_THROW_EXCEPTION(LogicException(oss.str()));
}



void
SVLocus::
mergeNode(
    const NodeIndexType fromIndex,
    const NodeIndexType toIndex)
{
    using namespace illumina::common;

    assert(fromIndex != toIndex);

#ifdef DEBUG_SVL
    static const std::string logtag("SVLocus::mergeNode");
    log_os << logtag << " from: " << fromIndex << " to: " << toIndex << " size: " << size() << "\n";
#endif

    SVLocusNode& fromNode(getNode(fromIndex));
    SVLocusNode& toNode(getNode(toIndex));

#ifdef DEBUG_SVL
    log_os << logtag << " BEFORE fromNode: " << fromNode;
    log_os << logtag << " BEFORE toNode: " << toNode;
#endif

    if (fromNode.interval.tid != toNode.interval.tid)
    {
        std::ostringstream oss;
        oss << "ERROR: Attempting to merge nodes on different chromosomes\n"
            << "\tNode1: " << fromNode
            << "\tNode2: " << toNode;
        BOOST_THROW_EXCEPTION(LogicException(oss.str()));
    }

    // store node relationship before merging regions:
    const bool isFromRegionRightmost(toNode.interval.range < fromNode.interval.range);

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

    // now take all fromNode edges and 'redirect' them to the toNode index
    //
    BOOST_FOREACH(edges_type::value_type& fromNodeEdgeIter, fromNode)
    {
        // alias value_type components (not required, but makes the logic easier to follow):
        const NodeIndexType& fromNodeEdgeIndex(fromNodeEdgeIter.first);
        SVLocusEdge& fromNodeEdge(fromNodeEdgeIter.second);

#ifdef DEBUG_SVL
        // is this edge between the to and from nodes?
        const bool isToFromEdge(fromNodeEdgeIndex == toIndex);

        log_os << logtag << " handle fromEdge: " << _index << ":" << fromNodeEdgeIndex << " isToFromEdge: " << isToFromEdge << "\n";
#endif


        // is this a self edge of the from node?
        const bool isSelfFromEdge(fromNodeEdgeIter.first == fromIndex);

        if (isSelfFromEdge)
        {
            // self-edge needs to be handled as a special case:
            toNode.mergeEdge(toIndex,fromNodeEdgeIter.second);
            continue;
        }

        // Check for the special case when there is an edge between from and to, in this case
        // the counts have to be handled so that counts in each region still approximate
        // fragment support. Normally (the chimera case) -- a single fragment will create
        // edges and nodes with weight one. If this is a non-chimera and the nodes collide and
        // merge, we want to prevent the evidence from being doubled when it should not be.
        //
        // To achieve this, we remove the edge counts from the 'right'-most region of the two nodes being merged. This is an
        // approximate solution, but very simple to add into the graph without blowing up per-node/edge storage.
        //
        const bool isFromToEdge(fromNodeEdgeIndex == toIndex);
        if (isFromToEdge)
        {
            SVLocusEdge* clearedEdgePtr(&fromNodeEdge);
            if (! isFromRegionRightmost)
            {
                edges_type::iterator toNodeEdgeIter(toNode.edges.find(fromIndex));
                assert(toNodeEdgeIter != toNode.edges.end());
                clearedEdgePtr=(&(toNodeEdgeIter->second));
            }
            toNode.count -= clearedEdgePtr->getCount();
            clearedEdgePtr->clearCount();
        }

        // update local edge:
        toNode.mergeEdge(fromNodeEdgeIndex,fromNodeEdge);

        // update remote inputNodeEdgeIter
        {
            SVLocusNode& remoteNode(getNode(fromNodeEdgeIndex));
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

            remoteNode.mergeEdge(toIndex, oldRemoteIter->second);
        }
    }

#ifdef DEBUG_SVL
    log_os << logtag << " AFTER toNode: " << toNode;
#endif

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
        if (edge.second.getCount() >= minMergeEdgeCount) return false;
        if (getEdge(edge.first,nodeIndex).getCount() >= minMergeEdgeCount) return false;
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
    static const std::string logtag("SVLocus::cleanNodeCore");
    log_os << logtag << " nodeAddy: " << _index << ":" << nodeIndex << "\n";
#endif

    unsigned totalCleaned(0);
    SVLocusNode& queryNode(getNode(nodeIndex));

    std::vector<NodeIndexType> eraseEdges;
    BOOST_FOREACH(edges_type::value_type& edgeIter, queryNode)
    {
        if (0 != edgeIter.second.getCount())
        {
            if (edgeIter.second.getCount() < minMergeEdgeCount)
            {
                // clean criteria met -- go ahead and erase edge count:
                assert(queryNode.count>=edgeIter.second.getCount());
                totalCleaned += edgeIter.second.getCount();
                queryNode.count -= edgeIter.second.getCount();
                edgeIter.second.clearCount();
            }
        }

        if (0 == edgeIter.second.getCount())
        {
            // if the out edge count is zero, see if the in-edge count is also zero --
            // if so, erase edge
            //
            const SVLocusEdge& fromRemoteEdge(getEdge(edgeIter.first,nodeIndex));
            if (0 == fromRemoteEdge.getCount())
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
#ifdef DEBUG_SVL
        log_os << logtag << " deleting edge: " << _index << ":" << nodeIndex << "->" << _index << ":" << toIndex << "\n";
#endif
        clearEdgePair(nodeIndex,toIndex);
    }

    // if true add the target node to the erase list:
    if ((0 == queryNode.edges.size()) && (0 == queryNode.count))
    {
        emptyNodes.insert(nodeIndex);
    }

#ifdef DEBUG_SVL
    log_os << logtag << " emptyEdges:\n";
    BOOST_FOREACH(const NodeIndexType toIndex, eraseEdges)
    {
        log_os << logtag << "\tedge: " << _index << ":" << nodeIndex << "->" << _index << ":" << toIndex << "\n";
    }

    log_os << "cleanNodeCore emptyNodes\n";
    BOOST_FOREACH(const NodeIndexType nodeIndex2, emptyNodes)
    {
        log_os << logtag << "\tnodeAddy: " << _index << ":" << nodeIndex2 << "\n";
    }

    log_os << logtag << " totalCleaned: " << totalCleaned << "\n";
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

#ifdef DEBUG_SVL
    static const std::string logtag("SVLocus::clearNodeEdges");
    log_os << logtag << " from nodeIndex: " << nodePtr << "\n";
#endif

    SVLocusNode& node(getNode(nodePtr));
    BOOST_FOREACH(edges_type::value_type& edgeIter, node)
    {

#ifdef DEBUG_SVL
        log_os << logtag << " clearing remote Index: " << edgeIter.first << "\n";
#endif
        // skip self edge (otherwise we invalidate iterators in this foreach loop)
        if (edgeIter.first == nodePtr) continue;

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

#ifdef DEBUG_SVL
        log_os << logtag << " remote clearing Index: " << thisRemoteIter->first << "\n";
#endif
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
    static const std::string logtag("SVLocus::eraseNode");
    log_os << logtag << " " << _index << ":" << nodePtr << " transfer_in: " << _index << ":" << fromPtr << " \n";

    log_os << logtag << " BEFORE: " << getNode(nodePtr) << "\n";
#endif

    if (fromPtr != nodePtr)
    {
#ifdef DEBUG_SVL
        log_os << logtag << " transfer_in: BEFORE: " << getNode(fromPtr) << "\n";
#endif
        // reassign fromNode's remote edges before shifting its address:
        //
        bool isHandleSelfEdge(false);
        SVLocusNode& fromNode(getNode(fromPtr));
        BOOST_FOREACH(edges_type::value_type& edgeIter, fromNode)
        {
            const bool isSelfEdge(edgeIter.first == fromPtr);

            if (isSelfEdge)
            {
                isHandleSelfEdge=true;
                continue;
            }

            SVLocusNode& remoteNode(getNode(edgeIter.first));
            edges_type& remoteEdges(remoteNode.edges);
            remoteEdges.insert(std::make_pair(nodePtr,getEdge(edgeIter.first,fromPtr)));
            remoteEdges.erase(fromPtr);
        }

        if (isHandleSelfEdge)
        {
            fromNode.edges.insert(std::make_pair(nodePtr,getEdge(fromPtr,fromPtr)));
            fromNode.edges.erase(fromPtr);
        }

        notifyDelete(nodePtr);
        _graph[nodePtr] = _graph[fromPtr];
        notifyAdd(nodePtr);

#ifdef DEBUG_SVL
        log_os << logtag << " transfer_in: AFTER: " << getNode(nodePtr) << "\n";
#endif
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
           << " out_count: " << edgeIter.second.getCount()
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
        sum += getEdge(edgeIter.first,nodeIndex).getCount();
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
           << " out_count: " << edgeIter.second.getCount()
           << " in_count: " << getEdge(edgeIter.first,nodeIndex).getCount() << "\n";
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
mergeSelfOverlap()
{
    const unsigned nodeSize(size());
    for (unsigned nodeIndex(0); nodeIndex<nodeSize; ++nodeIndex)
    {
        for (unsigned nodeIndex2(nodeIndex+1); nodeIndex2<nodeSize; ++nodeIndex2)
        {
            const unsigned revNodeIndex(nodeSize-(nodeIndex+1));
            const unsigned revNodeIndex2(nodeSize-(nodeIndex2+1));
            SVLocusNode& node1(getNode(revNodeIndex));
            SVLocusNode& node2(getNode(revNodeIndex2));

            // test whether 1 and 2 intersect, if they do, merge this into a self-edge node:
            if (! node2.interval.isIntersect(node1.interval)) continue;
            mergeNode(revNodeIndex,revNodeIndex2);
            eraseNode(revNodeIndex);
            break;
        }
    }
}



void
SVLocus::
checkState(const bool isCheckConnected) const
{
    using namespace illumina::common;

    const unsigned nodeSize(size());
    if (0 == nodeSize) return;

    for (unsigned nodeIndex(0); nodeIndex<nodeSize; ++nodeIndex)
    {
        const SVLocusNode& node(getNode(nodeIndex));

        // check that that every edge has a return path:
        BOOST_FOREACH(const edges_type::value_type& edgeIter, node)
        {
            getEdge(edgeIter.first,nodeIndex);
        }

        // check that node and edge counts are consistent:
        unsigned edgeCount(0);
        BOOST_FOREACH(const edges_type::value_type& edgeIter, node)
        {
            edgeCount += edgeIter.second.getCount();
        }

        if (edgeCount != node.count)
        {
            std::ostringstream oss;
            oss << "ERROR: SVLocusNode " << _index << ":" << nodeIndex << " has inconsistent counts. NodeCount: " << node.count << " EdgeCount: " << edgeCount << "\n";
            oss << "\tnode: " << node;
            BOOST_THROW_EXCEPTION(LogicException(oss.str()));
        }
    }

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
