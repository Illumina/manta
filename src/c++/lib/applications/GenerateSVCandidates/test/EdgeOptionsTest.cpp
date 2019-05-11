//
// Manta - Structural Variant and Indel Caller
// Copyright (c) 2013-2019 Illumina, Inc.
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

#include "boost/test/unit_test.hpp"

#include "EdgeOptions.hpp"

BOOST_AUTO_TEST_SUITE(EdgeOptionsTest_test_suite)

// Test the defalut locus edge Options for SVLocusGraph
// edge iteration and noise edge filtration
BOOST_AUTO_TEST_CASE(test_LocusEdgeOptions)
{
  struct LocusEdgeOptions lo;
  BOOST_REQUIRE_EQUAL(lo.isNodeIndex1, false);
  BOOST_REQUIRE_EQUAL(lo.isNodeIndex2, false);
  BOOST_REQUIRE_EQUAL(lo.locusIndex, 0);
  BOOST_REQUIRE_EQUAL(lo.nodeIndex1, 0);
  BOOST_REQUIRE_EQUAL(lo.nodeIndex2, 0);
}

// Test the default edge options for SVLocusGraph edge
// iteration and noise edge filtration
BOOST_AUTO_TEST_CASE(test_EdgeOptions)
{
  struct EdgeOptions eo;
  BOOST_REQUIRE_EQUAL(eo.isLocusIndex, false);
  BOOST_REQUIRE_EQUAL(eo.binCount, 1);
  BOOST_REQUIRE_EQUAL(eo.binIndex, 0);
  BOOST_REQUIRE_EQUAL(eo.graphNodeMaxEdgeCount, 10);
}

BOOST_AUTO_TEST_SUITE_END()
