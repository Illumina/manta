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

#pragma once

#include "blt_util/observer.hh"
#include "blt_util/known_pos_range2.hh"

#include "boost/foreach.hpp"
#include "boost/serialization/map.hpp"
#include "boost/serialization/vector.hpp"
#include "boost/serialization/split_member.hpp"

#include <iostream>
#include <limits>
#include <map>
#include <vector>


//#define DEBUG_SVL


#ifdef DEBUG_SVL
#include "blt_util/log.hh"

#include <iostream>
#endif


struct SVLocusSet;



// all internal locations use a chromosome index number
struct GenomeInterval
{
    GenomeInterval(
        const int32_t initTid = 0,
        const pos_t beginPos = 0,
        const pos_t endPos = 0) :
        tid(initTid),
        range(beginPos,endPos)
    {}

    /// does this intersect a second GenomeInterval?
    bool
    isIntersect(const GenomeInterval& gi) const
    {
        if (tid != gi.tid) return false;
        return range.is_range_intersect(gi.range);
    }

    bool
    operator<(const GenomeInterval& rhs) const
    {
        if (tid<rhs.tid) return true;
        if (tid == rhs.tid)
        {
            return (range<rhs.range);
        }
        return false;
    }

    bool
    operator==(const GenomeInterval& rhs) const
    {
        return ((tid==rhs.tid) && (range==rhs.range));
    }

    template<class Archive>
    void serialize(Archive& ar, const unsigned /* version */)
    {
        ar& tid& range;
    }

    int32_t tid;
    known_pos_range2 range;
};


std::ostream&
operator<<(std::ostream& os, const GenomeInterval& gi);

BOOST_CLASS_IMPLEMENTATION(GenomeInterval, boost::serialization::object_serializable)



struct SVLocusEdge
{
    SVLocusEdge(const unsigned init_count = 0) :
        count(init_count)
    {}

    // merge edge into this one
    //
    void
    mergeEdge(SVLocusEdge& edge)
    {
        count += edge.count;
    }

    template<class Archive>
    void serialize(Archive& ar, const unsigned /* version */)
    {
        ar& count;
    }

    unsigned short count;
};


std::ostream&
operator<<(std::ostream& os, const SVLocusEdge& edge);

BOOST_CLASS_IMPLEMENTATION(SVLocusEdge, boost::serialization::object_serializable)



typedef unsigned NodeIndexType;



struct SVLocusNode
{
    typedef std::map<NodeIndexType,SVLocusEdge> edges_type;
    typedef edges_type::iterator iterator;
    typedef edges_type::const_iterator const_iterator;

    SVLocusNode() :
        count(0)
    {}

    // specialized copy ctor which offsets all address:
    SVLocusNode(
        const SVLocusNode& in,
        const unsigned offset) :
        count(in.count),
        interval(in.interval)
    {
        BOOST_FOREACH(const edges_type::value_type& val, in)
        {
            edges.insert(std::make_pair(val.first+offset,val.second));
        }
    }

    unsigned
    size() const
    {
        return edges.size();
    }

    iterator
    begin()
    {
        return edges.begin();
    }

    iterator
    end()
    {
        return edges.end();
    }

    const_iterator
    begin() const
    {
        return edges.begin();
    }

    const_iterator
    end() const
    {
        return edges.end();
    }

    unsigned
    outCount() const
    {
        unsigned sum(0);
        BOOST_FOREACH(const edges_type::value_type& edgeIter, *this)
        {
            sum += edgeIter.second.count;
        }
        return sum;
    }

    template<class Archive>
    void serialize(Archive& ar,const unsigned /* version */)
    {
        ar& count& interval& edges;
    }

    friend std::ostream&
    operator<<(std::ostream& os, const SVLocusNode& node);

    //////////////////  data:
    unsigned short count;
    GenomeInterval interval;

    edges_type edges;
};


std::ostream&
operator<<(std::ostream& os, const SVLocusNode& node);

BOOST_CLASS_IMPLEMENTATION(SVLocusNode, boost::serialization::object_serializable)



typedef unsigned LocusIndexType;
typedef std::pair<bool, std::pair<LocusIndexType,NodeIndexType> > SVLocusNodeMoveMessage;




/// \brief a set of regions containing dependent SV evidence
///
/// An SV locus is a region hypothetically containing the breakends of 1 to many
/// SVs.
///
/// The locus is composed of a set of non-overlapping contiguous genomic regions and links between them.
/// Each link has an associated evidence count.
///
/// This class internally manages the node shared pointers in a synced data structure, there's probably a better
/// way to do this with transform_iterator, but I've always regretted using that.
///
struct SVLocus : public notifier<SVLocusNodeMoveMessage>
{
    typedef std::vector<SVLocusNode> graph_type;

    typedef graph_type::iterator iterator;
    typedef graph_type::const_iterator const_iterator;

    typedef SVLocusNode::edges_type edges_type;

    friend struct SVLocusSet;


    SVLocus() :
        _index(0)
    {}

    bool
    empty() const
    {
        return _graph.empty();
    }

    unsigned
    size() const
    {
        return _graph.size();
    }

    iterator
    begin()
    {
        return _graph.begin();
    }

    iterator
    end()
    {
        return _graph.end();
    }

    const_iterator
    begin() const
    {
        return _graph.begin();
    }

    const_iterator
    end() const
    {
        return _graph.end();
    }

    LocusIndexType
    getIndex() const
    {
        return _index;
    }

    const SVLocusNode&
    getNode(const NodeIndexType nodePtr) const
    {
        assert(nodePtr<_graph.size());
        return _graph[nodePtr];
    }

    NodeIndexType
    addNode(
        const GenomeInterval interval,
        const unsigned count = 1)
    {
        NodeIndexType nodePtr(newGraphNode());
        SVLocusNode& node(getNode(nodePtr));
        node.interval = interval;
        node.count=count;
        notifyAdd(nodePtr);
        return nodePtr;
    }

    NodeIndexType
    addRemoteNode(
        const GenomeInterval interval)
    {
        return addNode(interval,0);
    }

    // an edge count is only added on on from->to
    void
    linkNodes(
        const NodeIndexType fromIndex,
        const NodeIndexType toIndex,
        const unsigned fromCount = 1,
        const unsigned toCount = 0)
    {
        SVLocusNode& fromNode(getNode(fromIndex));
        SVLocusNode& toNode(getNode(toIndex));
        assert(0 == fromNode.edges.count(toIndex));
        assert(0 == toNode.edges.count(fromIndex));
        fromNode.edges.insert(std::make_pair(toIndex,SVLocusEdge(fromCount)));
        toNode.edges.insert(std::make_pair(fromIndex,SVLocusEdge(toCount)));
    }

    /// find all node indices connected to startIndex
    ///
    /// non-recursive version
    void
    findConnected(
            const NodeIndexType startIndex,
            std::set<NodeIndexType>& connected) const;

    /// debug func to check that internal data-structures are in
    /// a consistent state
    void
    checkState(const bool isCheckConnected = false) const;

    // total the evidence count of all in-edges to this node
    unsigned
    getNodeInCount(const LocusIndexType nodeIndex) const;

    // a fancier version of the SVLocusNode dumper which can
    // report in-edge information
    void
    dumpNode(
            std::ostream& os,
            const LocusIndexType nodeIndex) const;

    template<class Archive>
    void save(Archive& ar, const unsigned /* version */) const
    {
        ar << _graph;
    }

    template<class Archive>
    void load(Archive& ar, const unsigned /* version */)
    {
        clear();
        ar >> _graph;
    }

    BOOST_SERIALIZATION_SPLIT_MEMBER()

private:

    SVLocusNode&
    getNode(const NodeIndexType nodePtr)
    {
        assert(nodePtr<_graph.size());
        return _graph[nodePtr];
    }

    void
    updateIndex(const LocusIndexType& index)
    {
        _index=index;
    }

    // return true if node contains no out or in edges greater than
    // minMergeEdgeCount
    bool
    isNoiseNode(
            const unsigned minMergeEdgeCount,
            const NodeIndexType nodeIndex) const;

    /// remove all unmerged noise in-edges of node and
    /// return list of nodes which should be deleted
    void
    cleanNodeCore(
            const unsigned minMergeEdgeCount,
            const NodeIndexType nodeIndex,
            std::set<NodeIndexType>& emptyNodes);

    /// remove all unmerged noise in-edges of node and possibly
    /// delete empty nodes
    void
    cleanNode(
            const unsigned minMergeEdgeCount,
            const NodeIndexType nodeIndex);

    /// remove all unmerged noise edges and nodes
    void
    clean(const unsigned minMergeEdgeCount);

    void
    clearNodeEdges(const NodeIndexType nodePtr);

    void
    getEdgeException(
            const NodeIndexType fromIndex,
            const NodeIndexType toIndex) const;

    /// return from->to edge
    const SVLocusEdge&
    getEdge(const NodeIndexType fromIndex,
            const NodeIndexType toIndex) const
    {
        const SVLocusNode& fromNode(getNode(fromIndex));
        edges_type::const_iterator i(fromNode.edges.find(toIndex));
        if(i == fromNode.edges.end()) getEdgeException(fromIndex,toIndex);
        return i->second;
    }

    /// erase edges in both directions:
    void
    clearEdgePair(
            const NodeIndexType index1,
            const NodeIndexType index2)
    {
        clearEdge(index1,index2);
        clearEdge(index2,index1);
    }

    /// erase edge in one direction
    void
    clearEdge(
            const NodeIndexType fromIndex,
            const NodeIndexType toIndex)
    {
        getNode(fromIndex).edges.erase(toIndex);
    }

    /// copy fromLocus into this locus (this should be an intermediate part of a locus merge)
    void
    copyLocus(const SVLocus& fromLocus)
    {
        assert(&fromLocus != this);

        // simple method -- copy everything in with an offset in all index numbers:
        const unsigned offset(_graph.size());
        BOOST_FOREACH(const SVLocusNode& fromNode, fromLocus)
        {
            const NodeIndexType nodeIndex(newGraphNode());
            getNode(nodeIndex) = SVLocusNode(fromNode, offset);
            notifyAdd(nodeIndex);
        }
    }

    /// join from node into to node
    ///
    /// from node is effectively destroyed,
    //// because all incoming edges will be updated
    ///
    void
    mergeNode(
        const NodeIndexType fromIndex,
        const NodeIndexType toIndex);

    // remove node
    //
    void
    eraseNode(const NodeIndexType nodePtr);

    // remove a list of node ids
    void
    eraseNodes(const std::set<NodeIndexType>& nodes);

    void
    clear()
    {
        for (NodeIndexType i(0); i<size(); ++i) notifyDelete(i);
        _graph.clear();
    }

    NodeIndexType
    newGraphNode()
    {
        static const unsigned maxIndex(std::numeric_limits<NodeIndexType>::max());
        unsigned index(_graph.size());
        assert(index<maxIndex);
        _graph.resize(index+1);
        return static_cast<NodeIndexType>(index);
    }

    void
    notifyAdd(const NodeIndexType nodePtr)
    {
#ifdef DEBUG_SVL
        log_os << "SVLocusNotifier: Add node: " << _index << ":" << nodePtr << "\n";
#endif
        notify_observers(std::make_pair(true,std::make_pair(_index,nodePtr)));
    }

    void
    notifyDelete(const NodeIndexType nodePtr)
    {
#ifdef DEBUG_SVL
        log_os << "SVLocusNotifier: Delete node: " << _index << ":" << nodePtr << "\n";
#endif
        notify_observers(std::make_pair(false,std::make_pair(_index,nodePtr)));
    }

    graph_type _graph;
    LocusIndexType _index;
};


std::ostream&
operator<<(std::ostream& os, const SVLocus& locus);

BOOST_CLASS_IMPLEMENTATION(SVLocus, boost::serialization::object_serializable)
