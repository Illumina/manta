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


BOOST_AUTO_TEST_CASE( test_SVLocus1 ) {

    // construct a simple two-node locus
    SVLocus locus1;
    SVLocusNode* nodePtr = locus1.addNode(1,1000,2000);
    locus1.addNode(1,3000,4000,nodePtr);

    BOOST_CHECK_EQUAL(locus1.size(),2u);

    BOOST_FOREACH(const SVLocusNode* nodePtr1, locus1)
    {
        BOOST_CHECK_EQUAL(nodePtr1->edgeSize(),1u);
    }
}


BOOST_AUTO_TEST_CASE( test_SVLocusNodeMerge) {
    SVLocus locus1;
    SVLocusNode* nodePtr1 = locus1.addNode(1,1000,2000);

    SVLocus locus2;
    SVLocusNode* nodePtr2 = locus2.addNode(1,1500,2500);

    nodePtr1->mergeNode(*nodePtr2);

    BOOST_CHECK_EQUAL(nodePtr1->count,2u);
    BOOST_CHECK_EQUAL(nodePtr1->interval.range.begin_pos,1000u);
    BOOST_CHECK_EQUAL(nodePtr1->interval.range.end_pos,2500u);
}


BOOST_AUTO_TEST_SUITE_END()

