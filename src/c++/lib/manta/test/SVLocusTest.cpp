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


#include "boost/foreach.hpp"
#include "boost/test/unit_test.hpp"

#include "manta/SVLocus.hh"


BOOST_AUTO_TEST_SUITE( test_SVLocus )

BOOST_AUTO_TEST_CASE( test_GenomeInterval ) {

    // test that GenomeInterval sorting follows expect:
    std::vector<GenomeInterval> test;

    test.push_back(GenomeInterval(1,15,19));
    test.push_back(GenomeInterval(1,15,22));
    test.push_back(GenomeInterval(1,10,20));
    test.push_back(GenomeInterval(2,5,10));
    test.push_back(GenomeInterval(2,8,10));

    std::sort(test.begin(),test.end());

    BOOST_REQUIRE_EQUAL(test[0],GenomeInterval(1,10,20));
    BOOST_REQUIRE_EQUAL(test[2],GenomeInterval(1,15,22));
    BOOST_REQUIRE_EQUAL(test[4],GenomeInterval(2,8,10));
}


BOOST_AUTO_TEST_CASE( test_SVLocus1 ) {

    // construct a simple two-node locus
    SVLocus locus1;
    SVLocusNode* nodePtr = locus1.addNode(1,10,20);
    locus1.addNode(1,30,40,nodePtr);

    BOOST_REQUIRE_EQUAL(locus1.size(),2u);

    BOOST_FOREACH(const SVLocusNode* nodePtr1, locus1)
    {
        BOOST_REQUIRE_EQUAL(nodePtr1->edgeSize(),1u);
    }
}


BOOST_AUTO_TEST_CASE( test_SVLocusNodeMerge) {
    SVLocus locus1;
    SVLocusNode* nodePtr1 = locus1.addNode(1,10,20);

    SVLocus locus2;
    SVLocusNode* nodePtr2 = locus2.addNode(1,15,25);

    nodePtr1->mergeNode(*nodePtr2);

    BOOST_REQUIRE_EQUAL(nodePtr1->count,2u);
    BOOST_REQUIRE_EQUAL(nodePtr1->interval.range.begin_pos,10);
    BOOST_REQUIRE_EQUAL(nodePtr1->interval.range.end_pos,25);
}


BOOST_AUTO_TEST_CASE( test_SVLocusClearEdges ) {

    // construct a diamond four-node locus
    //
    //  1
    // 2 3
    //  4
    //
    SVLocus locus1;
    SVLocusNode* nodePtr1 = locus1.addNode(1,10,20);
    SVLocusNode* nodePtr2 = locus1.addNode(1,30,40,nodePtr1);
    SVLocusNode* nodePtr3 = locus1.addNode(1,50,60,nodePtr1);
    SVLocusNode* nodePtr4 = locus1.addNode(1,70,80,nodePtr2);
    nodePtr4->addEdge(*nodePtr3);

    // now disconnect 1 from 2,3:
    nodePtr1->clearEdges();

    BOOST_REQUIRE_EQUAL(locus1.size(),4u);

    BOOST_REQUIRE_EQUAL(nodePtr1->edgeSize(),0u);
    BOOST_REQUIRE_EQUAL(nodePtr2->edgeSize(),1u);
    BOOST_REQUIRE_EQUAL(nodePtr3->edgeSize(),1u);
    BOOST_REQUIRE_EQUAL(nodePtr4->edgeSize(),2u);
}


BOOST_AUTO_TEST_SUITE_END()

