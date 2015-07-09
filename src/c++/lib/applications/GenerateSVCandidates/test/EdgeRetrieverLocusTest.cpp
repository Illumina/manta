// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Manta - Structural Variant and Indel Caller
// Copyright (c) 2013-2015 Illumina, Inc.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//

///
/// \author Chris Saunders
///

#include "boost/test/unit_test.hpp"

#include "applications/GenerateSVCandidates/EdgeRetrieverLocus.hh"
#include "svgraph/SVLocusSet.hh"

#include "svgraph/test/SVLocusTestUtil.hh"


BOOST_AUTO_TEST_SUITE( test_EdgeRetrieverLocus )


BOOST_AUTO_TEST_CASE( test_EdgeRetrieverLocusSimple )
{
    SVLocus locus1;
    locusAddPair(locus1,1,10,20,2,30,40);

    SVLocus locus2;
    locusAddPair(locus2,3,10,20,4,30,40);

    SVLocusSetOptions sopt;
    sopt.minMergeEdgeObservations = 1;
    SVLocusSet set1(sopt);
    set1.merge(locus1);
    set1.merge(locus2);
    set1.checkState(true,true);

    BOOST_REQUIRE_EQUAL(set1.size(), 2u);

    LocusEdgeOptions lopt;
    lopt.locusIndex = 0;
    EdgeRetrieverLocus edger(set1, 0, lopt);

    BOOST_REQUIRE( edger.next() );

    EdgeInfo edge = edger.getEdge();
    BOOST_REQUIRE_EQUAL(edge.locusIndex, 0u);
    BOOST_REQUIRE_EQUAL(edge.nodeIndex1, 0u);
    BOOST_REQUIRE_EQUAL(edge.nodeIndex2, 1u);

    BOOST_REQUIRE(! edger.next() );
}


BOOST_AUTO_TEST_SUITE_END()

