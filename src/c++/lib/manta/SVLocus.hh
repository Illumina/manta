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

#include "blt_util/pos_range.hh"

#include "boost/foreach.hpp"
#include "boost/serialization/map.hpp"
#include "boost/serialization/set.hpp"
#include <boost/serialization/split_member.hpp>
#include "boost/shared_ptr.hpp"

#include <iostream>
#include <map>
#include <set>
#include <vector>


//#define DEBUG_SVL



// all internal locations use a chromosome index number
// and zero indexed position:
struct GenomeLocation
{
    GenomeLocation() :
        tid(0),
        pos(0)
    {}

    int32_t tid;
    int32_t pos;
};



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
    void serialize(Archive & ar, const unsigned /* version */)
    {
        ar & tid & range;
    }

    int32_t tid;
    known_pos_range range;
};


std::ostream&
operator<<(std::ostream& os, const GenomeInterval& gi);

BOOST_CLASS_IMPLEMENTATION(GenomeInterval, boost::serialization::object_serializable)



struct SVLocusNode;


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
    void serialize(Archive & ar, const unsigned /* version */)
    {
        ar & count;
    }

    unsigned count;
};


std::ostream&
operator<<(std::ostream& os, const SVLocusEdge& edge);

BOOST_CLASS_IMPLEMENTATION(SVLocusEdge, boost::serialization::object_serializable)



struct SVLocusNode
{
    typedef std::map<SVLocusNode*,SVLocusEdge> edges_type;
    typedef edges_type::const_iterator const_iterator;

    SVLocusNode() :
        count(0)
    {}

    unsigned
    size() const
    {
        return _edges.size();
    }

    const_iterator
    begin() const
    {
        return _edges.begin();
    }

    const_iterator
    end() const
    {
        return _edges.end();
    }

    void
    addEdge(SVLocusNode& linkTo,
            const bool isMakeReciprical = true)
    {
        // no self edges allowed:
        if (&linkTo == this) return;

        _edges.insert(std::make_pair(&linkTo,SVLocusEdge(1)));

        if (! isMakeReciprical) return;
        linkTo.addEdge(*this,false);
    }

    /// join the argument node into this one
    ///
    /// inputNode is required to intersect this node. inputNode is effectively destroyed,
    //// because all incoming edges will be updated
    ///
    void
    mergeNode(SVLocusNode& node);

    /// remove all outgoing and incoming edges to this node
    void
    clearEdges();

    /// debug function to check consistency
    void
    checkState() const;


    template<class Archive>
    void serialize(Archive & ar,const unsigned /* version */)
    {
        ar & count & interval & _edges;
    }

    friend std::ostream&
    operator<<(std::ostream& os, const SVLocusNode& node);

    //////////////////  data:
    unsigned count;
    GenomeInterval interval;

private:
    edges_type _edges;
};


std::ostream&
operator<<(std::ostream& os, const SVLocusNode& node);

BOOST_CLASS_IMPLEMENTATION(SVLocusNode, boost::serialization::object_serializable)



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
struct SVLocus
{
    typedef std::set<SVLocusNode*> graph_type;

    typedef graph_type::iterator iterator;
    typedef graph_type::const_iterator const_iterator;


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

    SVLocusNode*
    addNode(
        const int32_t tid,
        const int32_t beginPos,
        const int32_t endPos,
        SVLocusNode* linkTo = NULL)
    {
        SVLocusNode* nodePtr(newGraphNode());
        nodePtr->interval.tid=tid;
        nodePtr->interval.range.set_range(beginPos,endPos);
        nodePtr->count+=1;

        if (NULL != linkTo) {
            nodePtr->addEdge(*linkTo);
        }
        return nodePtr;
    }

    // copy a node from fromLocus into this locus
    void
    copyNode(const SVLocus& fromLocus,
             SVLocusNode* nodePtr)
    {
        assert(&fromLocus != this);

        const shared_map& fromSmap(fromLocus._smap);
        shared_map::const_iterator fromIter(fromSmap.find(nodePtr));
        assert(fromIter != fromSmap.end());

        shared_map::const_iterator toIter(_smap.find(nodePtr));
        assert(toIter == _smap.end());

        _graph.insert(nodePtr);
        _smap.insert(std::make_pair(nodePtr,fromIter->second));
    }

    // remove node
    //
    void
    erase(SVLocusNode* nodePtr)
    {
        assert(NULL != nodePtr);
        iterator iter(_graph.find(nodePtr));
        if (iter == _graph.end()) return;

        shared_map::iterator siter(_smap.find(nodePtr));
        assert(siter != _smap.end());

        nodePtr->clearEdges();

        _graph.erase(iter);
        _smap.erase(siter);
    }

    void
    clear()
    {
        _graph.clear();
        _smap.clear();
    }

    /// debug func to check that internal data-structures are in
    /// a consistent state
    void
    checkState() const;

    template<class Archive>
    void save(Archive & ar, const unsigned /* version */) const
    {
        // boost::serialize has shared_ptr options, but I'd rather not figure it out,
        // instead we will recreate shared_map on load:
        ar << _graph;
    }

    template<class Archive>
    void load(Archive & ar, const unsigned /* version */)
    {
        // boost::serialize has shared_ptr options, but I'd rather not figure it out,
        // instead we will recreate shared_map on load:
        clear();
        ar >> _graph;
        BOOST_FOREACH(SVLocusNode* nodePtr, *this)
        {
            _smap[nodePtr] = boost::shared_ptr<SVLocusNode>(nodePtr);
        }
    }

    BOOST_SERIALIZATION_SPLIT_MEMBER()

private:

    SVLocusNode*
    newGraphNode()
    {
        SVLocusNode* nodePtr(new SVLocusNode());
        boost::shared_ptr<SVLocusNode> sPtr(nodePtr);
        _graph.insert(nodePtr);
        _smap.insert(std::make_pair(nodePtr,sPtr));
        return nodePtr;
    }

    typedef boost::shared_ptr<SVLocusNode> shared_type;
    typedef std::map<const SVLocusNode*,shared_type> shared_map;

    graph_type _graph;
    shared_map _smap;
};


std::ostream&
operator<<(std::ostream& os, const SVLocus& locus);

BOOST_CLASS_IMPLEMENTATION(SVLocus, boost::serialization::object_serializable)
