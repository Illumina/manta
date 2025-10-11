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

#include "svgraph/SVLocus.hpp"
#include "test/testSVLocusUtil.hpp"

BOOST_AUTO_TEST_SUITE(SVLocusNode_test_suite)

BOOST_AUTO_TEST_CASE(test_SVLocusNode_EdgeManager)
{
  // test the new edge manager for SVLocusNode
  SVLocus       locus1;
  NodeIndexType nodePtr1     = locus1.addNode(GenomeInterval(1, 10, 20));
  NodeIndexType nodePtr1copy = locus1.addNode(GenomeInterval(1, 10, 20));
  locus1.linkNodes(nodePtr1, nodePtr1copy, 1, 1);
  locus1.mergeSelfOverlap();

  BOOST_REQUIRE_EQUAL(locus1.size(), 1u);

  const SVLocusNode& node(static_cast<const SVLocus&>(locus1).getNode(0));

  const SVLocusEdgeManager em = node.getEdgeManager();

  BOOST_REQUIRE_EQUAL(em.getMap().size(), 1u);

  BOOST_REQUIRE_EQUAL(em.getMap().begin()->first, 0u);
  BOOST_REQUIRE_EQUAL(em.getMap().begin()->second.getCount(), 1u);
}

BOOST_AUTO_TEST_SUITE_END()
