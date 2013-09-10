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

#include "applications/GenerateSVCandidates/EdgeRetrieverBin.hh"
#include "svgraph/SVLocusSet.hh"

#include "svgraph/test/SVLocusTestUtil.hh"

#include <iostream>


BOOST_AUTO_TEST_SUITE( test_EdgeRetrieverBin )


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

    EdgeRetrieverBin edger(set1, 0, 1, 0);

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
    const NodeIndexType nodePtr1 = locus1.addNode(GenomeInterval(1,10,20),3);
    const NodeIndexType nodePtr2 = locus1.addRemoteNode(GenomeInterval(2,30,40));
    locus1.linkNodes(nodePtr1,nodePtr2);
    const NodeIndexType nodePtr3 = locus1.addRemoteNode(GenomeInterval(3,30,40));
    locus1.linkNodes(nodePtr1,nodePtr3);
    const NodeIndexType nodePtr4 = locus1.addRemoteNode(GenomeInterval(4,30,40));
    locus1.linkNodes(nodePtr1,nodePtr4);
    SVLocus locus4;
    locusAddPair(locus4,7,10,20,8,30,40);

    SVLocusSet set1(1);
    set1.merge(locus1);
    set1.merge(locus4);
    set1.checkState(true,true);

    static const unsigned binTotal(2);
    for (unsigned binIndex(0); binIndex<binTotal; ++binIndex)
    {
        EdgeRetrieverBin edger(set1, 0, binTotal, binIndex);

        BOOST_REQUIRE( edger.next() );

        EdgeInfo edge = edger.getEdge();

        if (binIndex == 0)
        {
            BOOST_REQUIRE_EQUAL(edge.locusIndex, 0u);
            BOOST_REQUIRE_EQUAL(edge.nodeIndex1, 0u);
            BOOST_REQUIRE_EQUAL(edge.nodeIndex2, 1u);
        }
        else
        {
            BOOST_REQUIRE_EQUAL(edge.locusIndex, 0u);
            BOOST_REQUIRE_EQUAL(edge.nodeIndex1, 0u);
            BOOST_REQUIRE_EQUAL(edge.nodeIndex2, 3u);
        }
        BOOST_REQUIRE( edger.next() );

        edge = edger.getEdge();
        if (binIndex == 0)
        {
            BOOST_REQUIRE_EQUAL(edge.locusIndex, 0u);
            BOOST_REQUIRE_EQUAL(edge.nodeIndex1, 0u);
            BOOST_REQUIRE_EQUAL(edge.nodeIndex2, 2u);
        }
        else
        {
            BOOST_REQUIRE_EQUAL(edge.locusIndex, 1u);
            BOOST_REQUIRE_EQUAL(edge.nodeIndex1, 0u);
            BOOST_REQUIRE_EQUAL(edge.nodeIndex2, 1u);
        }

        BOOST_REQUIRE( ! edger.next() );
    }
}


BOOST_AUTO_TEST_CASE( test_EdgeRetrieverManyBin2 )
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
        EdgeRetrieverBin edger(set1, 0, 3, binIndex);

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
        EdgeRetrieverBin edger(set1, 0, 3, binIndex);

        while (edger.next())
        {
            count++;
        }
    }

    BOOST_REQUIRE_EQUAL(count,7u);
}



BOOST_AUTO_TEST_CASE( test_EdgeRetrieverOddBinSelfEdge )
{
    SVLocus locus1;
    locusAddPair(locus1,1,10,20,1,10,20);
    SVLocus locus2;
    locusAddPair(locus2,3,10,20,3,10,20);
    SVLocus locus3;
    locusAddPair(locus3,5,10,20,5,10,20);
    SVLocus locus4;
    locusAddPair(locus4,7,10,20,8,30,40);
    SVLocus locus5;
    locusAddPair(locus5,9,10,20,9,10,20);
    SVLocus locus6;
    locusAddPair(locus6,11,10,20,12,30,40);
    SVLocus locus7;
    locusAddPair(locus7,13,10,20,14,30,40);

    locus1.mergeSelfOverlap();
    locus2.mergeSelfOverlap();
    locus3.mergeSelfOverlap();
    locus5.mergeSelfOverlap();

    SVLocusSet set1(1);
    set1.merge(locus1);
    set1.merge(locus2);
    set1.merge(locus3);
    set1.merge(locus4);
    set1.merge(locus5);
    set1.merge(locus6);
    set1.merge(locus7);
    set1.checkState(true,true);

    std::set<unsigned> loci;
    unsigned count(0);
    for (unsigned binIndex(0); binIndex<3; ++binIndex)
    {
        EdgeRetrieverBin edger(set1, 0, 3, binIndex);

        while (edger.next())
        {
            const EdgeInfo& edge(edger.getEdge());
            BOOST_REQUIRE_EQUAL(loci.count(edge.locusIndex),0u);
            loci.insert(edge.locusIndex);
            count++;
        }
    }

    BOOST_REQUIRE_EQUAL(count,7u);
}


BOOST_AUTO_TEST_SUITE_END()

