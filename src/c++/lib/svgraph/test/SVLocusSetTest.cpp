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

#include "svgraph/SVLocusSet.hh"

#include "SVLocusTestUtil.hh"


BOOST_AUTO_TEST_SUITE( test_SVLocusSet )


BOOST_AUTO_TEST_CASE( test_SVLocusMerge )
{
    // test merge of overlapping loci

    // construct a simple two-node locus
    SVLocus locus1;
    locusAddPair(locus1,1,10,20,2,30,40);

    SVLocus locus2;
    locusAddPair(locus2,1,10,20,2,30,40);

    SVLocusSet set1(2);
    set1.merge(locus1);
    set1.merge(locus2);
    set1.checkState(true,true);
    const SVLocusSet& cset1(set1);

    BOOST_REQUIRE_EQUAL(cset1.nonEmptySize(),1u);
    BOOST_REQUIRE_EQUAL(cset1.getLocus(0).size(),2u);
}



BOOST_AUTO_TEST_CASE( test_SVLocusMultiOverlapMerge )
{
    // test merge of overlapping loci, reproduces production failure

    SVLocus locus1;
    locusAddPair(locus1,1,10,20,12,30,40);

    SVLocus locus2;
    locusAddPair(locus2,2,10,20,12,50,60);

    SVLocus locus3;
    locusAddPair(locus3,3,10,20,12,35,55);


    SVLocusSet set1(1);
    set1.merge(locus1);
    set1.merge(locus2);
    set1.merge(locus3);
    set1.checkState(true,true);
    const SVLocusSet& cset1(set1);

    GenomeInterval testInterval(12,30,60);

    BOOST_REQUIRE_EQUAL(cset1.nonEmptySize(),1u);
    BOOST_REQUIRE_EQUAL(cset1.getLocus(0).size(),4u);

    bool isFound(false);
    BOOST_FOREACH(const SVLocusNode& node, cset1.getLocus(0))
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
        NodeIndexType nodePtr1 = locus1.addNode(GenomeInterval(1,10,20));
        NodeIndexType nodePtr2 = locus1.addRemoteNode(GenomeInterval(1,30,40));
        NodeIndexType nodePtr3 = locus1.addRemoteNode(GenomeInterval(1,50,60));
        locus1.linkNodes(nodePtr1,nodePtr2);
        locus1.linkNodes(nodePtr1,nodePtr3);
    }

    SVLocus locus2;
    locusAddPair(locus2,1,10,60,2,10,60);

    SVLocusSet set1(1);
    set1.merge(locus1);
    set1.merge(locus2);
    set1.checkState(true,true);
    const SVLocusSet& cset1(set1);

    GenomeInterval testInterval(1,10,60);

    BOOST_REQUIRE_EQUAL(cset1.nonEmptySize(),1u);
    BOOST_REQUIRE_EQUAL(cset1.getLocus(0).size(),2u);

    bool isFound(false);
    BOOST_FOREACH(const SVLocusNode& node, cset1.getLocus(0))
    {
        if (node.interval == testInterval) isFound=true;
    }
    BOOST_REQUIRE(isFound);
}


BOOST_AUTO_TEST_CASE( test_SVLocusMultiOverlapMerge3 )
{
    // test merge of overlapping loci, reproduces production failure

    SVLocus locus1;
    locusAddPair(locus1,1,10,20,3,10,20);

    SVLocus locus2;
    locusAddPair(locus2,1,30,40,4,10,20);

    SVLocus locus3;
    locusAddPair(locus3,2,30,40,5,10,20);

    SVLocus locus4;
    locusAddPair(locus4,1,15,35,6,10,20);

    SVLocus locus5;
    locusAddPair(locus5,2,15,35,7,10,20);

    SVLocusSet set1(1);
    set1.merge(locus1);
    set1.merge(locus2);
    set1.merge(locus3);
    set1.merge(locus4);
    set1.merge(locus5);
    set1.checkState(true,true);
    const SVLocusSet& cset1(set1);

    GenomeInterval testInterval(1,10,40);

    BOOST_REQUIRE_EQUAL(cset1.nonEmptySize(),2u);
    BOOST_REQUIRE_EQUAL(cset1.getLocus(0).size(),4u);

    bool isFound(false);
    BOOST_FOREACH(const SVLocusNode& node, cset1.getLocus(0))
    {
        if (node.interval == testInterval) isFound=true;
    }
    BOOST_REQUIRE(isFound);
}



BOOST_AUTO_TEST_CASE( test_SVLocusMultiOverlapMerge4 )
{
    // test merge of overlapping loci, reproduces production failure

    SVLocus locus1;
    locusAddPair(locus1,1,10,60,2,20,30);

    SVLocus locus2;
    locusAddPair(locus2,1,40,50,1,20,30);

    SVLocusSet set1(1);
    set1.merge(locus1);
    set1.merge(locus2);

    const SVLocusSet& cset1(set1);
    cset1.checkState(true,true);


    GenomeInterval testInterval(1,10,60);

    BOOST_REQUIRE_EQUAL(cset1.nonEmptySize(),1u);
    BOOST_REQUIRE_EQUAL(cset1.getLocus(0).size(),2u);

    bool isFound(false);
    BOOST_FOREACH(const SVLocusNode& node, cset1.getLocus(0))
    {
        if (node.interval == testInterval) isFound=true;
    }
    BOOST_REQUIRE(isFound);
}



BOOST_AUTO_TEST_CASE( test_SVLocusNoiseMerge )
{
    SVLocus locus1;
    locusAddPair(locus1,1,10,60,2,20,30);

    SVLocus locus2;
    locusAddPair(locus2,1,10,60,2,20,30);

    SVLocus locus3;
    locusAddPair(locus3,1,10,60,3,20,30);

    {
        SVLocusSet set1(1);
        set1.merge(locus1);
        set1.merge(locus2);
        set1.merge(locus3);
        const SVLocusSet& cset1(set1);

        cset1.checkState(true,true);
        BOOST_REQUIRE_EQUAL(cset1.nonEmptySize(),1u);
        BOOST_REQUIRE_EQUAL(cset1.getLocus(0).size(),3u);
    }

    {
        SVLocusSet set1(2);
        set1.merge(locus1);
        set1.merge(locus2);
        set1.merge(locus3);
        const SVLocusSet& cset1(set1);

        cset1.checkState(true,true);
        BOOST_REQUIRE_EQUAL(cset1.nonEmptySize(),2u);
        BOOST_REQUIRE_EQUAL(cset1.getLocus(0).size(),2u);
    }

    {
        SVLocusSet set1(3);
        set1.merge(locus1);
        set1.merge(locus2);
        set1.merge(locus3);
        const SVLocusSet& cset1(set1);

        cset1.checkState(true,true);
        BOOST_REQUIRE_EQUAL(cset1.nonEmptySize(),3u);
        BOOST_REQUIRE_EQUAL(cset1.getLocus(0).size(),2u);
    }
}



BOOST_AUTO_TEST_CASE( test_SVLocusNoiseClean )
{
    SVLocus locus1;
    locusAddPair(locus1,1,10,60,2,20,30);

    SVLocus locus2;
    locusAddPair(locus2,1,10,60,2,20,30);

    SVLocus locus3;
    locusAddPair(locus3,1,10,60,3,20,30);

    {
        SVLocusSet set1(2);
        set1.merge(locus1);
        set1.merge(locus2);
        set1.merge(locus3);
        const SVLocusSet& cset1(set1);

        cset1.checkState(true,true);
        BOOST_REQUIRE_EQUAL(cset1.nonEmptySize(),2u);
        BOOST_REQUIRE_EQUAL(cset1.getLocus(0).size(),2u);

        set1.clean();

        BOOST_REQUIRE_EQUAL(cset1.nonEmptySize(),1u);
        BOOST_REQUIRE_EQUAL(cset1.getLocus(0).size(),2u);
    }

    {
        SVLocusSet set1(2);
        set1.merge(locus1);
        set1.merge(locus2);
        set1.merge(locus3);
        const SVLocusSet& cset1(set1);

        BOOST_REQUIRE_EQUAL(cset1.nonEmptySize(),2u);
        BOOST_REQUIRE_EQUAL(cset1.getLocus(1).size(),2u);

        set1.cleanRegion(GenomeInterval(3,0,70));

        BOOST_REQUIRE_EQUAL(cset1.nonEmptySize(),2u);
        BOOST_REQUIRE_EQUAL(cset1.getLocus(1).size(),2u);

        set1.cleanRegion(GenomeInterval(1,0,70));

        BOOST_REQUIRE_EQUAL(cset1.nonEmptySize(),1u);
        BOOST_REQUIRE_EQUAL(cset1.getLocus(0).size(),2u);
    }

}



BOOST_AUTO_TEST_CASE( test_SVLocusNoiseCleanOrder )
{
    SVLocus locus1;
    locusAddPair(locus1,1,10,60,2,20,30);

    SVLocus locus2;
    locusAddPair(locus2,1,10,60,2,20,30);

    SVLocus locus3;
    locusAddPair(locus3,1,10,60,3,20,30);

    SVLocus locus4;
    locusAddPair(locus4,1,10,60,4,20,30);

    SVLocus locus5;
    locusAddPair(locus5,1,10,60,4,20,30);

    SVLocus locus6;
    locusAddPair(locus6,1,10,60,5,20,30);

    {
        SVLocusSet set1(2);
        set1.merge(locus1);
        set1.merge(locus2);
        set1.merge(locus3);
        set1.merge(locus4);
        set1.merge(locus5);
        set1.merge(locus6);
        const SVLocusSet& cset1(set1);

        BOOST_REQUIRE_EQUAL(cset1.nonEmptySize(),3u);
        BOOST_REQUIRE_EQUAL(cset1.getLocus(1).size(),2u);

        set1.cleanRegion(GenomeInterval(1,0,70));

        BOOST_REQUIRE_EQUAL(cset1.nonEmptySize(),1u);
        BOOST_REQUIRE_EQUAL(cset1.getLocus(0).size(),3u);
    }

}



BOOST_AUTO_TEST_CASE( test_SVLocusNoiseCleanRemote )
{
    SVLocus locus1;
    locusAddPair(locus1,1,100,110,1,10,20);

    {
        SVLocusSet set1(2);
        set1.merge(locus1);
        const SVLocusSet& cset1(set1);

        BOOST_REQUIRE_EQUAL(cset1.nonEmptySize(),1u);
        BOOST_REQUIRE_EQUAL(cset1.getLocus(0).size(),2u);

        set1.cleanRegion(GenomeInterval(1,0,120));

        BOOST_REQUIRE_EQUAL(cset1.nonEmptySize(),0u);
    }

}



BOOST_AUTO_TEST_CASE( test_SVLocusEvidenceRange )
{
    SVLocus locus1;
    {
        NodeIndexType node1 = locus1.addNode(GenomeInterval(1,100,110));
        NodeIndexType node2 = locus1.addRemoteNode(GenomeInterval(2,100,110));
        locus1.linkNodes(node1,node2);
        locus1.setNodeEvidence(node1,known_pos_range2(50,60));
    }

    SVLocus locus2;
    {
        NodeIndexType node1 = locus2.addNode(GenomeInterval(1,100,110));
        NodeIndexType node2 = locus2.addRemoteNode(GenomeInterval(2,100,110));
        locus2.linkNodes(node1,node2);
        locus2.setNodeEvidence(node1,known_pos_range2(30,40));
    }

    {
        SVLocusSet set1(2);
        set1.merge(locus1);
        set1.merge(locus2);
        const SVLocusSet& cset1(set1);

        BOOST_REQUIRE_EQUAL(cset1.nonEmptySize(),1u);
        BOOST_REQUIRE_EQUAL(cset1.getLocus(0).getNode(0).evidenceRange,known_pos_range2(30,60));
    }

}






BOOST_AUTO_TEST_CASE( test_SVLocusNoiseOverlap )
{
    // adapted from a production failure case:
    SVLocus locus1;
    locusAddPair(locus1,1,10,60,2,20,30);
    SVLocus locus2;
    locusAddPair(locus2,1,10,60,2,20,30);
    SVLocus locus3;
    locusAddPair(locus3,1,59,70,3,20,30);
    SVLocus locus4;
    locusAddPair(locus4,1,65,70,3,20,30);

    {
        SVLocusSet set1(2);
        set1.merge(locus1);
        set1.merge(locus2);
        set1.merge(locus3);
        set1.merge(locus4);
        const SVLocusSet& cset1(set1);

        set1.finalize();
        cset1.checkState(true,true);
    }
}


BOOST_AUTO_TEST_SUITE_END()

