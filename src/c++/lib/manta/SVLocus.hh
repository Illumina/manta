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
    SVLocusEdge() :
        count(0),
        next(NULL)
    {}

    unsigned count;
    SVLocusNode* next;
};



struct SVLocusNode
{

    SVLocusNode() :
        count(0)
    {}

    unsigned count;
    GenomeInterval interval;
    std::vector<SVLocusEdge> edges;
};



/// \brief a set of regions containing dependent SV evidence
///
/// An SV locus is a region hypothetically containing the breakends of 1 to many
/// SVs.
///
/// The locus is composed of a set of non-overlapping contiguous genomic regions and links between them.
/// Each link has an associated evidence count.
///
///
struct SVLocus {

    bool
    empty() const
    {
        return graph.empty();
    }


    SVLocusNode*
    addNode(
            const int32_t tid,
            const int32_t begin,
            const int32_t end,
            SVLocusNode* linkTo = NULL)
    {
        SVLocusNode* nodePtr(new SVLocusNode());
        boost::shared_ptr<SVLocusNode> sPtr(nodePtr);
        graph.push_back(sPtr);
        SVLocusNode& node (*(graph.back()));
        node.interval.tid=tid;
        node.interval.range.set_range(begin,end);
        node.count+=1;
        if(NULL != linkTo) {
            node.edges.push_back(SVLocusEdge());
            node.edges.back().count=1;
            node.edges.back().next=linkTo;

            linkTo->edges.push_back(SVLocusEdge());
            linkTo->edges.back().count=1;
            linkTo->edges.back().next=nodePtr;
        }
        return nodePtr;
    }

    void
    clear()
    {
        graph.clear();
    }


    std::vector<boost::shared_ptr<SVLocusNode> > graph;
};

