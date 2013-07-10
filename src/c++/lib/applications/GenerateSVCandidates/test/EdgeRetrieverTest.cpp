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

#include "applications/GenerateSVCandidates/EdgeRetriever.hh"
#include "svgraph/SVLocusSet.hh"

#include "svgraph/test/SVLocusTestUtil.hh"


BOOST_AUTO_TEST_SUITE( test_EdgeRetriever )


BOOST_AUTO_TEST_CASE( test_EdgeRetrieverOneBin )
{
    SVLocus locus1;
    locusAddPair(locus1,1,10,20,2,30,40);

    SVLocus locus2;
    locusAddPair(locus2,3,10,20,4,30,40);

    SVLocusSet set1(1);
    set1.merge(locus1);
    set1.merge(locus2);
    set1.checkState(true,true);

    EdgeRetriever edger(set1,1,0);

    BOOST_REQUIRE( edger.next() );

    EdgeInfo edge = edger.getEdge();
    BOOST_REQUIRE_EQUAL(edge.locusIndex, 0u);
    BOOST_REQUIRE_EQUAL(edge.nodeIndex1, 0u);
    BOOST_REQUIRE_EQUAL(edge.nodeIndex2, 1u);

    BOOST_REQUIRE( edger.next() );

    edge = edger.getEdge();
    BOOST_REQUIRE_EQUAL(edge.locusIndex, 1u);
    BOOST_REQUIRE_EQUAL(edge.nodeIndex1, 0u);
    BOOST_REQUIRE_EQUAL(edge.nodeIndex2, 1u);

    BOOST_REQUIRE( ! edger.next() );
}



BOOST_AUTO_TEST_CASE( test_EdgeRetrieverManyBin )
{
    SVLocus locus1;
    locusAddPair(locus1,1,10,20,2,30,40);
    SVLocus locus2;
    locusAddPair(locus2,3,10,20,4,30,40);
    SVLocus locus3;
    locusAddPair(locus3,5,10,20,6,30,40);
    SVLocus locus4;
    locusAddPair(locus4,7,10,20,8,30,40);
    SVLocus locus5;
    locusAddPair(locus5,9,10,20,10,30,40);
    SVLocus locus6;
    locusAddPair(locus6,11,10,20,12,30,40);

    SVLocusSet set1(1);
    set1.merge(locus1);
    set1.merge(locus2);
    set1.merge(locus3);
    set1.merge(locus4);
    set1.merge(locus5);
    set1.merge(locus6);
    set1.checkState(true,true);

    for (unsigned binIndex(0); binIndex<3; ++binIndex)
    {
        EdgeRetriever edger(set1,3,binIndex);

        BOOST_REQUIRE( edger.next() );

        EdgeInfo edge = edger.getEdge();
        BOOST_REQUIRE_EQUAL(edge.locusIndex, 0u+(binIndex*2));
        BOOST_REQUIRE_EQUAL(edge.nodeIndex1, 0u);
        BOOST_REQUIRE_EQUAL(edge.nodeIndex2, 1u);

        BOOST_REQUIRE( edger.next() );

        edge = edger.getEdge();
        BOOST_REQUIRE_EQUAL(edge.locusIndex, 1u+(binIndex*2));
        BOOST_REQUIRE_EQUAL(edge.nodeIndex1, 0u);
        BOOST_REQUIRE_EQUAL(edge.nodeIndex2, 1u);

        BOOST_REQUIRE( ! edger.next() );
    }
}



BOOST_AUTO_TEST_CASE( test_EdgeRetrieverOddBin )
{
    SVLocus locus1;
    locusAddPair(locus1,1,10,20,2,30,40);
    SVLocus locus2;
    locusAddPair(locus2,3,10,20,4,30,40);
    SVLocus locus3;
    locusAddPair(locus3,5,10,20,6,30,40);
    SVLocus locus4;
    locusAddPair(locus4,7,10,20,8,30,40);
    SVLocus locus5;
    locusAddPair(locus5,9,10,20,10,30,40);
    SVLocus locus6;
    locusAddPair(locus6,11,10,20,12,30,40);
    SVLocus locus7;
    locusAddPair(locus7,13,10,20,14,30,40);

    SVLocusSet set1(1);
    set1.merge(locus1);
    set1.merge(locus2);
    set1.merge(locus3);
    set1.merge(locus4);
    set1.merge(locus5);
    set1.merge(locus6);
    set1.merge(locus7);
    set1.checkState(true,true);

    unsigned count(0);
    for (unsigned binIndex(0); binIndex<3; ++binIndex)
    {
        EdgeRetriever edger(set1,3,binIndex);

        while (edger.next())
        {
            count++;
        }
    }

    BOOST_REQUIRE_EQUAL(count,7u);
}


BOOST_AUTO_TEST_SUITE_END()

