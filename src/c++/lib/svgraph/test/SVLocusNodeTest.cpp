// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Manta
// Copyright (c) 2013-2014 Illumina, Inc.
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

#define private public
#include "svgraph/SVLocus.hh"

#include "SVLocusTestUtil.hh"


BOOST_AUTO_TEST_SUITE( test_SVLocusNode )


BOOST_AUTO_TEST_CASE( test_SVLocusNode_EM )
{
    // test the new edge manager for SVLocusNode
    SVLocus locus1;
    NodeIndexType nodePtr1 = locus1.addNode(GenomeInterval(1,10,20));
    NodeIndexType nodePtr1copy = locus1.addNode(GenomeInterval(1,10,20));
    locus1.linkNodes(nodePtr1,nodePtr1copy,1,1);
    locus1.mergeSelfOverlap();

    BOOST_REQUIRE_EQUAL(locus1.size(),1u);

    SVLocusNode& node(locus1.getNode(0));

    SVLocusEdgeManager em = node.getEdgeManager();

    BOOST_REQUIRE_EQUAL(em.getMap().size(),1u);

    BOOST_REQUIRE_EQUAL(em.getMap().begin()->first,0u);
    BOOST_REQUIRE_EQUAL(em.getMap().begin()->second.getCount(),1u);
}


BOOST_AUTO_TEST_SUITE_END()
