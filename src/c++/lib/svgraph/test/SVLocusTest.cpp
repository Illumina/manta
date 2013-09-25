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

#include "boost/foreach.hpp"
#include "boost/test/unit_test.hpp"

#define private public
#include "svgraph/SVLocus.hh"

#include "SVLocusTestUtil.hh"


BOOST_AUTO_TEST_SUITE( test_SVLocus )


BOOST_AUTO_TEST_CASE( test_SVLocus1 )
{

    // construct a simple two-node locus
    SVLocus locus1;
    locusAddPair(locus1,1,10,20,1,30,40);

    BOOST_REQUIRE_EQUAL(locus1.size(),2u);

    BOOST_FOREACH(const SVLocusNode& node, static_cast<const SVLocus&>(locus1))
    {
        BOOST_REQUIRE_EQUAL(node.size(),1u);
    }
}


BOOST_AUTO_TEST_CASE( test_SVLocusNodeMerge2)
{
    SVLocus locus1;
    NodeIndexType nodePtr1 = locus1.addNode(GenomeInterval(1,10,20));
    NodeIndexType nodePtr2 = locus1.addNode(GenomeInterval(1,15,25));
    NodeIndexType nodePtr3 = locus1.addNode(GenomeInterval(3,10,20));
    NodeIndexType nodePtr4 = locus1.addNode(GenomeInterval(4,10,20));
    locus1.linkNodes(nodePtr1,nodePtr3);
    locus1.linkNodes(nodePtr2,nodePtr4);

    locus1.mergeNode(nodePtr2,nodePtr1);

    const SVLocusNode& node1(locus1.getNode(nodePtr1));

    BOOST_REQUIRE_EQUAL(node1.getCount(),2u);
    BOOST_REQUIRE_EQUAL(node1.interval.range.begin_pos(),10);
    BOOST_REQUIRE_EQUAL(node1.interval.range.end_pos(),25);
    BOOST_REQUIRE_EQUAL(node1.size(),2u);
}


BOOST_AUTO_TEST_CASE( test_SVLocusNodeMergeSelfEdge)
{
    SVLocus locus1;
    NodeIndexType nodePtr1 = locus1.addNode(GenomeInterval(1,10,20));
    NodeIndexType nodePtr2 = locus1.addNode(GenomeInterval(1,15,25));
    locus1.linkNodes(nodePtr1,nodePtr2);
    locus1.mergeSelfOverlap();

    const SVLocusNode& node1(locus1.getNode(0));

    BOOST_REQUIRE_EQUAL(node1.getCount(),1u);
    BOOST_REQUIRE_EQUAL(node1.interval.range.begin_pos(),10);
    BOOST_REQUIRE_EQUAL(node1.interval.range.end_pos(),25);
    BOOST_REQUIRE_EQUAL(node1.size(),1u);

    // test that the single edge of the merged node is to self:
    BOOST_REQUIRE_EQUAL(node1.begin()->first,0u);
    BOOST_REQUIRE_EQUAL(node1.outCount(),1u);
}


BOOST_AUTO_TEST_CASE( test_SVLocusNodeMergeSelfEdgeReverse)
{
    SVLocus locus1;
    NodeIndexType nodePtr1 = locus1.addNode(GenomeInterval(1,10,20));
    NodeIndexType nodePtr2 = locus1.addNode(GenomeInterval(1,15,25));
    locus1.linkNodes(nodePtr2,nodePtr1);
    locus1.mergeSelfOverlap();

    const SVLocusNode& node1(locus1.getNode(0));

    BOOST_REQUIRE_EQUAL(node1.getCount(),0u);
    BOOST_REQUIRE_EQUAL(node1.interval.range.begin_pos(),10);
    BOOST_REQUIRE_EQUAL(node1.interval.range.end_pos(),25);
    BOOST_REQUIRE_EQUAL(node1.size(),1u);

    // test that the single edge of the merged node is to self:
    BOOST_REQUIRE_EQUAL(node1.begin()->first,0u);
    BOOST_REQUIRE_EQUAL(node1.outCount(),0u);
}


BOOST_AUTO_TEST_CASE( test_SVLocusNodeMergeMultiSelfEdge )
{
    SVLocus locus1;
    NodeIndexType nodePtr1 = locus1.addNode(GenomeInterval(1,10,20));
    NodeIndexType nodePtr1copy = locus1.addNode(GenomeInterval(1,10,20));
    NodeIndexType nodePtr2 = locus1.addNode(GenomeInterval(1,15,25));
    locus1.linkNodes(nodePtr1,nodePtr1copy,1,1);
    locus1.linkNodes(nodePtr1,nodePtr2);

    locus1.mergeSelfOverlap();

    const SVLocusNode& node1(locus1.getNode(0));

    BOOST_REQUIRE_EQUAL(node1.getCount(),2u);
    BOOST_REQUIRE_EQUAL(node1.interval.range.begin_pos(),10);
    BOOST_REQUIRE_EQUAL(node1.interval.range.end_pos(),25);
    BOOST_REQUIRE_EQUAL(node1.size(),1u);

    // test that the single edge of the merged node is to self:
    BOOST_REQUIRE_EQUAL(node1.begin()->first,0u);
    BOOST_REQUIRE_EQUAL(node1.outCount(),2u);
}



BOOST_AUTO_TEST_CASE( test_SVLocusClearEdges )
{
    // construct a diamond four-node locus
    //
    //  1
    // 2 3
    //  4
    //
    SVLocus locus1;
    NodeIndexType nodePtr1 = locus1.addNode(GenomeInterval(1,10,20));
    NodeIndexType nodePtr2 = locus1.addNode(GenomeInterval(1,30,40));
    NodeIndexType nodePtr3 = locus1.addNode(GenomeInterval(1,50,60));
    NodeIndexType nodePtr4 = locus1.addNode(GenomeInterval(1,70,80));
    locus1.linkNodes(nodePtr1,nodePtr2);
    locus1.linkNodes(nodePtr1,nodePtr3);
    locus1.linkNodes(nodePtr2,nodePtr4);
    locus1.linkNodes(nodePtr3,nodePtr4);

    // now disconnect 1 from 2,3:
    locus1.clearNodeEdges(nodePtr1);

    BOOST_REQUIRE_EQUAL(locus1.size(),4u);

    BOOST_REQUIRE_EQUAL(locus1.getNode(nodePtr1).size(),0u);
    BOOST_REQUIRE_EQUAL(locus1.getNode(nodePtr2).size(),1u);
    BOOST_REQUIRE_EQUAL(locus1.getNode(nodePtr3).size(),1u);
    BOOST_REQUIRE_EQUAL(locus1.getNode(nodePtr4).size(),2u);
}


BOOST_AUTO_TEST_SUITE_END()

