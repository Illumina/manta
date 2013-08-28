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

#include "applications/GenerateSVCandidates/EdgeRetrieverLocus.hh"
#include "svgraph/SVLocusSet.hh"

#include "svgraph/test/SVLocusTestUtil.hh"

#include <iostream>


BOOST_AUTO_TEST_SUITE( test_EdgeRetrieverLocus )


BOOST_AUTO_TEST_CASE( test_EdgeRetrieverLocusSimple )
{
    SVLocus locus1;
    locusAddPair(locus1,1,10,20,2,30,40);

    SVLocus locus2;
    locusAddPair(locus2,3,10,20,4,30,40);

    SVLocusSet set1(1);
    set1.merge(locus1);
    set1.merge(locus2);
    set1.checkState(true,true);

    BOOST_REQUIRE_EQUAL(set1.size(), 2u);

    EdgeRetrieverLocus edger(set1, 0, 0);

    BOOST_REQUIRE( edger.next() );

    EdgeInfo edge = edger.getEdge();
    BOOST_REQUIRE_EQUAL(edge.locusIndex, 0u);
    BOOST_REQUIRE_EQUAL(edge.nodeIndex1, 0u);
    BOOST_REQUIRE_EQUAL(edge.nodeIndex2, 1u);

    BOOST_REQUIRE(! edger.next() );
}


BOOST_AUTO_TEST_SUITE_END()

