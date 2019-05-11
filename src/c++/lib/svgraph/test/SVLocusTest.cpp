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

#include "svgraph/SVLocus.hpp"
#include "test/testSVLocusUtil.hpp"

BOOST_AUTO_TEST_SUITE(SVLocus_test_suite)

BOOST_AUTO_TEST_CASE(test_SVLocus1)
{
  // construct a simple two-node locus
  SVLocus locus1;
  locusAddPair(locus1, 1, 10, 20, 1, 30, 40);

  BOOST_REQUIRE_EQUAL(locus1.size(), 2u);

  for (const SVLocusNode& node : static_cast<const SVLocus&>(locus1)) {
    BOOST_REQUIRE_EQUAL(node.size(), 1u);
  }
}

BOOST_AUTO_TEST_CASE(test_SVLocusNodeMerge2)
{
  SVLocus       locus1;
  NodeIndexType nodePtr1 = locus1.addNode(GenomeInterval(1, 10, 20));
  NodeIndexType nodePtr2 = locus1.addNode(GenomeInterval(1, 15, 25));
  NodeIndexType nodePtr3 = locus1.addNode(GenomeInterval(3, 10, 20));
  NodeIndexType nodePtr4 = locus1.addNode(GenomeInterval(4, 10, 20));
  locus1.linkNodes(nodePtr1, nodePtr3);
  locus1.linkNodes(nodePtr2, nodePtr4);

  //locus1.mergeNode(nodePtr2, nodePtr1, nullptr);
  locus1.mergeSelfOverlap();

  const SVLocusNode& node1(static_cast<const SVLocus&>(locus1).getNode(nodePtr1));

  BOOST_REQUIRE_EQUAL(node1.outCount(), 2u);
  BOOST_REQUIRE_EQUAL(node1.getInterval().range.begin_pos(), 10);
  BOOST_REQUIRE_EQUAL(node1.getInterval().range.end_pos(), 25);
  BOOST_REQUIRE_EQUAL(node1.size(), 2u);
}

BOOST_AUTO_TEST_CASE(test_SVLocusNodeMergeSelfEdge)
{
  SVLocus       locus1;
  NodeIndexType nodePtr1 = locus1.addNode(GenomeInterval(1, 10, 20));
  NodeIndexType nodePtr2 = locus1.addNode(GenomeInterval(1, 15, 25));
  locus1.linkNodes(nodePtr1, nodePtr2);
  locus1.mergeSelfOverlap();

  const SVLocusNode& node1(static_cast<const SVLocus&>(locus1).getNode(0));

  BOOST_REQUIRE_EQUAL(node1.outCount(), 1u);
  BOOST_REQUIRE_EQUAL(node1.getInterval().range.begin_pos(), 10);
  BOOST_REQUIRE_EQUAL(node1.getInterval().range.end_pos(), 25);
  BOOST_REQUIRE_EQUAL(node1.size(), 1u);

  // test that the single edge of the merged node is to self:
  BOOST_REQUIRE_EQUAL(node1.getEdgeManager().getMap().begin()->first, 0u);
}

BOOST_AUTO_TEST_CASE(test_SVLocusNodeMergeSelfEdgeReverse)
{
  SVLocus       locus1;
  NodeIndexType nodePtr1 = locus1.addNode(GenomeInterval(1, 10, 20));
  NodeIndexType nodePtr2 = locus1.addNode(GenomeInterval(1, 15, 25));
  locus1.linkNodes(nodePtr2, nodePtr1);
  locus1.mergeSelfOverlap();

  const SVLocusNode& node1(static_cast<const SVLocus&>(locus1).getNode(0));

  BOOST_REQUIRE_EQUAL(node1.outCount(), 1u);
  BOOST_REQUIRE_EQUAL(node1.getInterval().range.begin_pos(), 10);
  BOOST_REQUIRE_EQUAL(node1.getInterval().range.end_pos(), 25);
  BOOST_REQUIRE_EQUAL(node1.size(), 1u);

  // test that the single edge of the merged node is to self:
  BOOST_REQUIRE_EQUAL(node1.getEdgeManager().getMap().begin()->first, 0u);
}

BOOST_AUTO_TEST_CASE(test_SVLocusNodeMergeMultiSelfEdge)
{
  SVLocus       locus1;
  NodeIndexType nodePtr1     = locus1.addNode(GenomeInterval(1, 10, 20));
  NodeIndexType nodePtr1copy = locus1.addNode(GenomeInterval(1, 10, 20));
  NodeIndexType nodePtr2     = locus1.addNode(GenomeInterval(1, 15, 25));
  locus1.linkNodes(nodePtr1, nodePtr1copy, 1, 1);
  locus1.linkNodes(nodePtr1, nodePtr2);

  locus1.mergeSelfOverlap();

  const SVLocusNode& node1(static_cast<const SVLocus&>(locus1).getNode(0));

  BOOST_REQUIRE_EQUAL(node1.outCount(), 2u);
  BOOST_REQUIRE_EQUAL(node1.getInterval().range.begin_pos(), 10);
  BOOST_REQUIRE_EQUAL(node1.getInterval().range.end_pos(), 25);
  BOOST_REQUIRE_EQUAL(node1.size(), 1u);

  // test that the single edge of the merged node is to self:
  BOOST_REQUIRE_EQUAL(node1.getEdgeManager().getMap().begin()->first, 0u);
}

BOOST_AUTO_TEST_CASE(test_SVLocusClearEdges)
{
  // construct a diamond four-node locus
  //
  //  1
  // 2 3
  //  4
  //
  SVLocus       locus1;
  NodeIndexType nodePtr1 = locus1.addNode(GenomeInterval(1, 10, 20));
  NodeIndexType nodePtr2 = locus1.addNode(GenomeInterval(1, 30, 40));
  NodeIndexType nodePtr3 = locus1.addNode(GenomeInterval(1, 50, 60));
  NodeIndexType nodePtr4 = locus1.addNode(GenomeInterval(1, 70, 80));
  locus1.linkNodes(nodePtr1, nodePtr2);
  locus1.linkNodes(nodePtr1, nodePtr3);
  locus1.linkNodes(nodePtr2, nodePtr4);
  locus1.linkNodes(nodePtr3, nodePtr4);

  // now disconnect 1 from 2,3:
  locus1.clearNodeEdges(nodePtr1);

  const SVLocus& clocus1(locus1);
  BOOST_REQUIRE_EQUAL(clocus1.size(), 4u);

  BOOST_REQUIRE_EQUAL(clocus1.getNode(nodePtr1).size(), 0u);
  BOOST_REQUIRE_EQUAL(clocus1.getNode(nodePtr2).size(), 1u);
  BOOST_REQUIRE_EQUAL(clocus1.getNode(nodePtr3).size(), 1u);
  BOOST_REQUIRE_EQUAL(clocus1.getNode(nodePtr4).size(), 2u);
}

BOOST_AUTO_TEST_SUITE_END()
