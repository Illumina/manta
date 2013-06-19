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

#include "blt_util/pos_range.hh"

#include "boost/shared_ptr.hpp"

#include <map>
#include <set>
#include <vector>


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
    GenomeInterval() :
        tid(0),
        range(0,0)
    {}

    bool
    operator<(const GenomeInterval& rhs) const
    {
        if(tid<rhs.tid) return true;
        if(tid==rhs.tid)
        {
            return (range<rhs.range);
        }
        return false;
    }

    int32_t tid;
    known_pos_range range;
};


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

    unsigned count;
};



struct SVLocusNode
{

    SVLocusNode() :
        count(0)
    {}

    void
    addEdge(SVLocusNode& linkTo,
            const bool isMakeReciprical = true)
    {
        // no self edges allowed:
        if(&linkTo == this) return;

        _edges.insert(std::make_pair(&linkTo,SVLocusEdge(1)));

        if(! isMakeReciprical) return;
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

    unsigned count;
    GenomeInterval interval;

    typedef std::map<SVLocusNode*,SVLocusEdge> edges_type;
    edges_type _edges;
};



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

        if(NULL != linkTo) {
            nodePtr->addEdge(*linkTo);
        }
        return nodePtr;
    }

    // copy a node from fromLocus into this locus
    void
    copyNode(const SVLocus& fromLocus,
             SVLocusNode* nodePtr)
    {
        const shared_map& from_smap(fromLocus._smap);
        shared_map::const_iterator iter(from_smap.find(nodePtr));
        assert(iter != from_smap.end());

        _graph.insert(nodePtr);
        _smap.insert(std::make_pair(nodePtr,iter->second));
    }

    // remove node
    //
    void
    erase(SVLocusNode* nodePtr)
    {
        assert(NULL != nodePtr);
        iterator iter(_graph.find(nodePtr));
        if(iter == _graph.end()) return;

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
    typedef std::map<graph_type::value_type,shared_type> shared_map;

    graph_type _graph;
    shared_map _smap;
};

