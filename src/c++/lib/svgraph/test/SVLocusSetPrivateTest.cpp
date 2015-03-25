// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Manta
// Copyright (c) 2013-2015 Illumina, Inc.
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

#include "boost/test/unit_test.hpp"

// hack to call private methods of SVLocusSet:
#pragma clang diagnostic ignored "-Wkeyword-macro"
#define private public

#include "svgraph/SVLocusSet.hh"

#include "SVLocusTestUtil.hh"


BOOST_AUTO_TEST_SUITE( test_SVLocusSetPrivate )


static
unsigned
testOverlap(
    SVLocusSet& locusSet,
    const int32_t tid,
    const int32_t beginPos,
    const int32_t endPos)
{
    std::set<SVLocusSet::NodeAddressType> intersect;
    locusSet.getRegionIntersect(GenomeInterval(tid,beginPos,endPos),intersect);
    return intersect.size();
}



BOOST_AUTO_TEST_CASE( test_SVLocusIntersect )
{
    // construct a simple two-node locus
    SVLocus locus1;
    locusAddPair(locus1,1,10,20,2,30,40);

    SVLocusSet set1;
    set1.merge(locus1);
    set1.checkState(true,true);

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


BOOST_AUTO_TEST_CASE( test_SVLocusCombine )
{
    // test reassigning the locus numbers of non-overlapping loci in a set:

    // construct a simple two-node locus
    SVLocus locus1;
    locusAddPair(locus1,1,10,20,2,30,40);

    SVLocus locus2;
    locusAddPair(locus2,3,10,20,4,30,40);

    SVLocus locus3;
    locusAddPair(locus3,5,10,20,6,30,40);

    SVLocusSetOptions sopt;
    sopt.minMergeEdgeObservations = 1;
    SVLocusSet set1(sopt);
    set1.merge(locus1);
    set1.merge(locus2);
    set1.merge(locus3);

    set1.checkState(true,true);

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

