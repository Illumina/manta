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

#include "EdgeRetrieverLocus.cpp"
#include "EdgeRetrieverLocus.hpp"

#include "options/SVLocusSetOptions.hpp"
#include "svgraph/SVLocus.hpp"
#include "test/testSVLocusUtil.hpp"

BOOST_AUTO_TEST_SUITE(EdgeRetrieverLocus_test_suite)

BOOST_AUTO_TEST_CASE(test_EdgeRetrieverLocusSimple)
{
  SVLocus locus1;
  locusAddPair(locus1, 1, 10, 20, 2, 30, 40);

  SVLocus locus2;
  locusAddPair(locus2, 3, 10, 20, 4, 30, 40);

  SVLocusSetOptions sopt;
  sopt.minMergeEdgeObservations = 1;
  SVLocusSet set1(sopt);
  set1.merge(locus1);
  set1.merge(locus2);
  set1.checkState(true, true);

  BOOST_REQUIRE_EQUAL(set1.size(), 2u);

  LocusEdgeOptions lopt;
  lopt.locusIndex = 0;
  EdgeRetrieverLocus edger(set1, 0, lopt);

  BOOST_REQUIRE(edger.next());
  EdgeInfo edge = edger.getEdge();
  BOOST_REQUIRE_EQUAL(edge.locusIndex, 0u);
  BOOST_REQUIRE_EQUAL(edge.nodeIndex1, 0u);
  BOOST_REQUIRE_EQUAL(edge.nodeIndex2, 1u);
  BOOST_REQUIRE(!edger.next());

  // when someone is interested to generate candidates for only the
  // edge from node1 to node2 in one locus.
  // In the following test case graph has 5 nodes:
  // Node-0: GenomeInterval:1 [10, 20)
  // Node-1: GenomeInterval:2 [30, 40)
  // Node-2: GenomeInterval:4 [30, 40)
  // Node-3: GenomeInterval:5 [10, 20)
  // Node-4: GenomeInterval:6 [60, 70)
  // Let's say we want only candidate in between Node-0 and Node-2
  SVLocus locus3;
  locusAddPair(locus3, 1, 10, 20, 2, 30, 40);
  SVLocus locus4;
  locusAddPair(locus4, 1, 10, 20, 4, 30, 40);
  SVLocus locus5;
  locusAddPair(locus5, 5, 10, 20, 4, 30, 40);
  SVLocus locus6;
  locusAddPair(locus6, 5, 10, 20, 6, 60, 70);
  SVLocusSet set2(sopt);
  set2.merge(locus3);
  set2.merge(locus4);
  set2.merge(locus5);
  set2.merge(locus6);
  set2.checkState(true, true);
  lopt.locusIndex   = 0;  // report this locus only
  lopt.isNodeIndex1 = true;
  lopt.isNodeIndex2 = true;
  // Required Node index
  lopt.nodeIndex1 = 0;
  lopt.nodeIndex2 = 2;
  EdgeRetrieverLocus edger1(set2, 4, lopt);

  BOOST_REQUIRE(edger1.next());
  edge = edger1.getEdge();
  BOOST_REQUIRE_EQUAL(edge.locusIndex, 0u);
  BOOST_REQUIRE_EQUAL(edge.nodeIndex1, 0u);
  BOOST_REQUIRE_EQUAL(edge.nodeIndex2, 2u);
  BOOST_REQUIRE(!edger1.next());

  // Candidates for all edges touching node1 (i.e. Node-0) are interested
  lopt.isNodeIndex1 = true;
  lopt.isNodeIndex2 = false;
  // Required Node index
  lopt.nodeIndex1 = 0;
  EdgeRetrieverLocus edger2(set2, 4, lopt);

  BOOST_REQUIRE(edger2.next());
  edge = edger2.getEdge();
  BOOST_REQUIRE_EQUAL(edge.locusIndex, 0u);
  // Node-0 ----> Node-1
  BOOST_REQUIRE_EQUAL(edge.nodeIndex1, 0u);
  BOOST_REQUIRE_EQUAL(edge.nodeIndex2, 1u);
  BOOST_REQUIRE(edger2.next());
  edge = edger2.getEdge();
  // Node-0 ----> Node-2
  BOOST_REQUIRE_EQUAL(edge.nodeIndex1, 0u);
  BOOST_REQUIRE_EQUAL(edge.nodeIndex2, 2u);
  BOOST_REQUIRE(!edger2.next());

  lopt.isNodeIndex2 = true;
  // Node-0 and Node-2 both node size = 2
  // But here max Node edge count = 1.
  // So it's a noisy node. It will not return any edge
  EdgeRetrieverLocus edger3(set2, 1, lopt);
  BOOST_REQUIRE(!edger3.next());
}

BOOST_AUTO_TEST_SUITE_END()
