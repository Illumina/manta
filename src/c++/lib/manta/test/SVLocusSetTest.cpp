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



unsigned
testOverlap(
    SVLocusSet& locusSet,
    const int32_t tid,
    const int32_t beginPos,
    const int32_t endPos)
{
    SVLocusSet::LocusSetIndexerType intersect(locusSet);
    locusSet.getRegionIntersect(tid,beginPos,endPos,intersect);
    return intersect.size();
}



BOOST_AUTO_TEST_CASE( test_SVLocusIntersect ) {

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



BOOST_AUTO_TEST_CASE( test_SVLocusMerge ) {

    // test merge of overlapping loci

    // construct a simple two-node locus
    SVLocus locus1;
    NodeIndexType nodePtr1 = locus1.addNode(1,10,20);
    NodeIndexType nodePtr2 = locus1.addNode(2,30,40);
    locus1.linkNodes(nodePtr1,nodePtr2);

    SVLocus locus2;
    NodeIndexType nodePtr3 = locus2.addNode(1,10,20);
    NodeIndexType nodePtr4 = locus2.addNode(2,30,40);
    locus2.linkNodes(nodePtr3,nodePtr4);

    SVLocusSet set1;
    set1.merge(locus1);
    set1.merge(locus2);

    BOOST_REQUIRE_EQUAL(set1.nonEmptySize(),1u);
    BOOST_REQUIRE_EQUAL(set1._loci[0].size(),2u);
}


BOOST_AUTO_TEST_CASE( test_SVLocusCombine ) {

    // test reassigning the locus numbers of non-overlapping loci in a set:

    // construct a simple two-node locus
    SVLocus locus1;
    NodeIndexType nodePtr1 = locus1.addNode(1,10,20);
    NodeIndexType nodePtr2 = locus1.addNode(2,30,40);
    locus1.linkNodes(nodePtr1,nodePtr2);

    SVLocus locus2;
    NodeIndexType nodePtr3 = locus2.addNode(3,10,20);
    NodeIndexType nodePtr4 = locus2.addNode(4,30,40);
    locus2.linkNodes(nodePtr3,nodePtr4);

    SVLocus locus3;
    NodeIndexType nodePtr5 = locus3.addNode(5,10,20);
    NodeIndexType nodePtr6 = locus3.addNode(6,30,40);
    locus3.linkNodes(nodePtr5,nodePtr6);

    SVLocusSet set1;
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

