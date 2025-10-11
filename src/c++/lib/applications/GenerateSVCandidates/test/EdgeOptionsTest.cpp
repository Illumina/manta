//
// Manta - Structural Variant and Indel Caller
// Copyright (c) 2013-2025 Illumina, Inc.
//
// This program is licensed under the terms of the Polyform strict license
//
// ***As far as the law allows, the software comes as is, without
// any warranty or condition, and the licensor will not be liable
// to you for any damages arising out of these terms or the use
// or nature of the software, under any kind of legal claim.***
//
// You should have received a copy of the PolyForm Strict License 1.0.0
// along with this program.  If not, see <https://polyformproject.org/licenses/strict/1.0.0>.
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
