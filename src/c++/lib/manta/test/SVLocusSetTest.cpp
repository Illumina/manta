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

#include "boost/test/unit_test.hpp"

// hack to call private methods of SVLocusSet:
#define private public
#include "manta/SVLocusSet.hh"


BOOST_AUTO_TEST_SUITE( test_SVLocusSet )

static
unsigned
testOverlap(
    SVLocusSet& locusSet,
    const int32_t tid,
    const int32_t beginPos,
    const int32_t endPos)
{
    std::set<SVLocusSet::NodeAddressType> intersect;
    locusSet.getRegionIntersect(tid,beginPos,endPos,intersect);
    return intersect.size();
}



BOOST_AUTO_TEST_CASE( test_SVLocusIntersect )
{
    // construct a simple two-node locus
    SVLocus locus1;
    NodeIndexType nodePtr1 = locus1.addNode(1,10,20);
    NodeIndexType nodePtr2 = locus1.addNode(2,30,40);
    locus1.linkNodes(nodePtr1,nodePtr2);

    SVLocusSet set1;
    set1.merge(locus1);

    // test for various intersections:

    // non-overlap test:
    BOOST_REQUIRE_EQUAL(testOverlap(set1,1,1,2),0u);

    // left-edge overlap:
    BOOST_REQUIRE_EQUAL(testOverlap(set1,1,9,11),1u);

    // right-edge overlap:
    BOOST_REQUIRE_EQUAL(testOverlap(set1,1,19,21),1u);

    // non-overlap:
    BOOST_REQUIRE_EQUAL(testOverlap(set1,1,29,31),0u);

    // non-overlap (diff tid):
    BOOST_REQUIRE_EQUAL(testOverlap(set1,2,9,11),0u);
}



BOOST_AUTO_TEST_CASE( test_SVLocusMerge )
{
    // test merge of overlapping loci

    // construct a simple two-node locus
    SVLocus locus1;
    NodeIndexType nodePtr1 = locus1.addNode(1,10,20);
    NodeIndexType nodePtr2 = locus1.addRemoteNode(2,30,40);
    locus1.linkNodes(nodePtr1,nodePtr2);

    SVLocus locus2;
    NodeIndexType nodePtr3 = locus2.addNode(1,10,20);
    NodeIndexType nodePtr4 = locus2.addRemoteNode(2,30,40);
    locus2.linkNodes(nodePtr3,nodePtr4);

    SVLocusSet set1(2);
    set1.merge(locus1);
    set1.merge(locus2);

    BOOST_REQUIRE_EQUAL(set1.nonEmptySize(),1u);
    BOOST_REQUIRE_EQUAL(set1.getLocus(0).size(),2u);
}



BOOST_AUTO_TEST_CASE( test_SVLocusMultiOverlapMerge )
{
    // test merge of overlapping loci, reproduces production failure

    SVLocus locus1;
    {
        NodeIndexType nodePtr1 = locus1.addNode(1,10,20);
        NodeIndexType nodePtr2 = locus1.addRemoteNode(12,30,40);
        locus1.linkNodes(nodePtr1,nodePtr2);
    }

    SVLocus locus2;
    {
        NodeIndexType nodePtr1 = locus2.addNode(2,10,20);
        NodeIndexType nodePtr2 = locus2.addRemoteNode(12,50,60);
        locus2.linkNodes(nodePtr1,nodePtr2);
    }

    SVLocus locus3;
    {
        NodeIndexType nodePtr1 = locus3.addNode(3,10,20);
        NodeIndexType nodePtr2 = locus3.addRemoteNode(12,35,55);
        locus3.linkNodes(nodePtr1,nodePtr2);
    }


    SVLocusSet set1(1);
    set1.merge(locus1);
    set1.merge(locus2);
    set1.merge(locus3);

    GenomeInterval testInterval(12,30,60);

    BOOST_REQUIRE_EQUAL(set1.nonEmptySize(),1u);
    BOOST_REQUIRE_EQUAL(set1.getLocus(0).size(),4u);

    bool isFound(false);
    BOOST_FOREACH(const SVLocusNode& node, set1.getLocus(0))
    {
        if (node.interval == testInterval) isFound=true;
    }
    BOOST_REQUIRE(isFound);
}



BOOST_AUTO_TEST_CASE( test_SVLocusMultiOverlapMerge2 )
{
    // test merge of overlapping loci, reproduces production failure

    SVLocus locus1;
    {
        NodeIndexType nodePtr1 = locus1.addNode(1,10,20);
        NodeIndexType nodePtr2 = locus1.addRemoteNode(1,30,40);
        NodeIndexType nodePtr3 = locus1.addRemoteNode(1,50,60);
        locus1.linkNodes(nodePtr1,nodePtr2);
        locus1.linkNodes(nodePtr1,nodePtr3);
    }

    SVLocus locus2;
    {
        NodeIndexType nodePtr1 = locus2.addNode(1,10,60);
        NodeIndexType nodePtr2 = locus2.addRemoteNode(2,10,60);
        locus2.linkNodes(nodePtr1,nodePtr2);
    }

    SVLocusSet set1(1);
    set1.merge(locus1);
    set1.merge(locus2);

    GenomeInterval testInterval(1,10,60);

    BOOST_REQUIRE_EQUAL(set1.nonEmptySize(),1u);
    BOOST_REQUIRE_EQUAL(set1.getLocus(0).size(),2u);

    bool isFound(false);
    BOOST_FOREACH(const SVLocusNode& node, set1.getLocus(0))
    {
        if (node.interval == testInterval) isFound=true;
    }
    BOOST_REQUIRE(isFound);
}



BOOST_AUTO_TEST_CASE( test_SVLocusMultiOverlapMerge3 )
{
    // test merge of overlapping loci, reproduces production failure

    SVLocus locus1;
    {
        NodeIndexType nodePtr1 = locus1.addNode(1,10,20);
        NodeIndexType nodePtr2 = locus1.addRemoteNode(3,10,20);
        locus1.linkNodes(nodePtr1,nodePtr2);
    }

    SVLocus locus2;
    {
        NodeIndexType nodePtr1 = locus2.addNode(1,30,40);
        NodeIndexType nodePtr2 = locus2.addRemoteNode(4,10,20);
        locus2.linkNodes(nodePtr1,nodePtr2);
    }

    SVLocus locus3;
    {
        NodeIndexType nodePtr1 = locus3.addNode(2,30,40);
        NodeIndexType nodePtr2 = locus3.addRemoteNode(5,10,20);
        locus3.linkNodes(nodePtr1,nodePtr2);
    }

    SVLocus locus4;
    {
        NodeIndexType nodePtr1 = locus4.addNode(1,15,35);
        NodeIndexType nodePtr2 = locus4.addRemoteNode(6,10,20);
        locus4.linkNodes(nodePtr1,nodePtr2);
    }

    SVLocus locus5;
    {
        NodeIndexType nodePtr1 = locus5.addNode(2,15,35);
        NodeIndexType nodePtr2 = locus5.addRemoteNode(7,10,20);
        locus5.linkNodes(nodePtr1,nodePtr2);
    }

    SVLocusSet set1(1);
    set1.merge(locus1);
    set1.merge(locus2);
    set1.merge(locus3);
    set1.merge(locus4);
    set1.merge(locus5);

    GenomeInterval testInterval(1,10,40);

    BOOST_REQUIRE_EQUAL(set1.nonEmptySize(),2u);
    BOOST_REQUIRE_EQUAL(set1.getLocus(0).size(),4u);

    bool isFound(false);
    BOOST_FOREACH(const SVLocusNode& node, set1.getLocus(0))
    {
        if (node.interval == testInterval) isFound=true;
    }
    BOOST_REQUIRE(isFound);
}



BOOST_AUTO_TEST_CASE( test_SVLocusMultiOverlapMerge4 )
{
    // test merge of overlapping loci, reproduces production failure

    SVLocus locus1;
    {
        NodeIndexType nodePtr1 = locus1.addNode(1,10,60);
        NodeIndexType nodePtr2 = locus1.addNode(2,20,30);
        locus1.linkNodes(nodePtr1,nodePtr2);
    }

    SVLocus locus2;
    {
        NodeIndexType nodePtr1 = locus2.addNode(1,40,50);
        NodeIndexType nodePtr2 = locus2.addNode(1,20,30);
        locus2.linkNodes(nodePtr1,nodePtr2);
    }

    SVLocusSet set1(1);
    set1.merge(locus1);
    set1.merge(locus2);

    set1.checkState(true);

    GenomeInterval testInterval(1,10,60);

    BOOST_REQUIRE_EQUAL(set1.nonEmptySize(),1u);
    BOOST_REQUIRE_EQUAL(set1.getLocus(0).size(),2u);

    bool isFound(false);
    BOOST_FOREACH(const SVLocusNode& node, set1.getLocus(0))
    {
        if (node.interval == testInterval) isFound=true;
    }
    BOOST_REQUIRE(isFound);
}



BOOST_AUTO_TEST_CASE( test_SVLocusCombine )
{
    // test reassigning the locus numbers of non-overlapping loci in a set:

    // construct a simple two-node locus
    SVLocus locus1;
    {
        NodeIndexType nodePtr1 = locus1.addNode(1,10,20);
        NodeIndexType nodePtr2 = locus1.addRemoteNode(2,30,40);
        locus1.linkNodes(nodePtr1,nodePtr2);
    }

    SVLocus locus2;
    {
        NodeIndexType nodePtr3 = locus2.addNode(3,10,20);
        NodeIndexType nodePtr4 = locus2.addRemoteNode(4,30,40);
        locus2.linkNodes(nodePtr3,nodePtr4);
    }

    SVLocus locus3;
    {
        NodeIndexType nodePtr5 = locus3.addNode(5,10,20);
        NodeIndexType nodePtr6 = locus3.addRemoteNode(6,30,40);
        locus3.linkNodes(nodePtr5,nodePtr6);
    }

    SVLocusSet set1(1);
    set1.merge(locus1);
    set1.merge(locus2);
    set1.merge(locus3);

    BOOST_REQUIRE_EQUAL(set1._loci[0].size(),2u);
    BOOST_REQUIRE_EQUAL(set1._loci[1].size(),2u);
    BOOST_REQUIRE_EQUAL(set1._loci[2].size(),2u);
    set1.combineLoci(0,0);
    set1.combineLoci(2,0);

    BOOST_REQUIRE_EQUAL(set1._loci[0].size(),4u);
    BOOST_REQUIRE_EQUAL(set1._loci[1].size(),2u);
    BOOST_REQUIRE_EQUAL(set1._loci[2].size(),0u);
}


BOOST_AUTO_TEST_SUITE_END()

