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
    SVLocusSet::ins_type intersect;

    // non-overlap test:
    SVLocus locus2;
    SVLocusNode* nodePtr1 = locus2.addNode(tid,beginPos,endPos);
    locusSet.getNodeIntersect(nodePtr1,intersect);
    return intersect.size();
}



BOOST_AUTO_TEST_CASE( test_SVLocusIntersect ) {

    // construct a simple two-node locus
    SVLocus locus1;
    SVLocusNode* nodePtr = locus1.addNode(1,10,20);
    locus1.addNode(2,30,40,nodePtr);

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

BOOST_AUTO_TEST_SUITE_END()

