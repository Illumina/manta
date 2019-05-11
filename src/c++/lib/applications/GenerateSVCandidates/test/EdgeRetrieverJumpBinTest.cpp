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

#include "EdgeRetrieverJumpBin.cpp"
#include "EdgeRetrieverJumpBin.hpp"

#include "svgraph/SVLocusSet.hpp"
#include "test/testSVLocusUtil.hpp"

BOOST_AUTO_TEST_SUITE(EdgeRetrieverJumpBin_test_suite)

BOOST_AUTO_TEST_CASE(test_EdgeRetrieverJumpOneBin)
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

  EdgeRetrieverJumpBin edger(set1, 0, 1, 0);

  BOOST_REQUIRE(edger.next());

  EdgeInfo edge = edger.getEdge();
  BOOST_REQUIRE_EQUAL(edge.locusIndex, 0u);
  BOOST_REQUIRE_EQUAL(edge.nodeIndex1, 0u);
  BOOST_REQUIRE_EQUAL(edge.nodeIndex2, 1u);

  BOOST_REQUIRE(edger.next());

  edge = edger.getEdge();
  BOOST_REQUIRE_EQUAL(edge.locusIndex, 1u);
  BOOST_REQUIRE_EQUAL(edge.nodeIndex1, 0u);
  BOOST_REQUIRE_EQUAL(edge.nodeIndex2, 1u);

  BOOST_REQUIRE(!edger.next());

  // Here max node edge threshold as 1.
  // But all the nodes have two edges.
  // So all the nodes are noisy. Api does not return
  // any node.
  SVLocusSet set2(sopt);
  SVLocus    locus3;
  locusAddPair(locus3, 1, 10, 20, 2, 30, 40);
  SVLocus locus4;
  locusAddPair(locus4, 1, 10, 20, 6, 30, 40);
  SVLocus locus5;
  locusAddPair(locus5, 7, 10, 20, 2, 30, 40);
  SVLocus locus6;
  locusAddPair(locus6, 7, 10, 20, 6, 30, 40);
  set2.merge(locus3);
  set2.merge(locus4);
  set2.merge(locus5);
  set2.merge(locus6);
  set2.checkState(true, true);
  // max node edge cutoff is 1.
  EdgeRetrieverJumpBin edger1(set2, 1, 1, 0);
  BOOST_REQUIRE(!edger1.next());
}

BOOST_AUTO_TEST_CASE(test_EdgeRetrieverJumpManyBin)
{
  SVLocus             locus1;
  const NodeIndexType nodePtr1 = locus1.addNode(GenomeInterval(1, 10, 20));
  const NodeIndexType nodePtr2 = locus1.addNode(GenomeInterval(2, 30, 40));
  locus1.linkNodes(nodePtr1, nodePtr2);
  const NodeIndexType nodePtr3 = locus1.addNode(GenomeInterval(3, 30, 40));
  locus1.linkNodes(nodePtr1, nodePtr3);
  const NodeIndexType nodePtr4 = locus1.addNode(GenomeInterval(4, 30, 40));
  locus1.linkNodes(nodePtr1, nodePtr4);
  SVLocus locus4;
  locusAddPair(locus4, 7, 10, 20, 8, 30, 40);

  SVLocusSetOptions sopt;
  sopt.minMergeEdgeObservations = 1;
  SVLocusSet set1(sopt);
  set1.merge(locus1);
  set1.merge(locus4);
  set1.checkState(true, true);

  static const unsigned binTotal(2);
  for (unsigned binIndex(0); binIndex < binTotal; ++binIndex) {
    EdgeRetrieverJumpBin edger(set1, 0, binTotal, binIndex);

    BOOST_REQUIRE(edger.next());

    EdgeInfo edge = edger.getEdge();

    if (binIndex == 0) {
      BOOST_REQUIRE_EQUAL(edge.locusIndex, 0u);
      BOOST_REQUIRE_EQUAL(edge.nodeIndex1, 0u);
      BOOST_REQUIRE_EQUAL(edge.nodeIndex2, 1u);
    } else {
      BOOST_REQUIRE_EQUAL(edge.locusIndex, 0u);
      BOOST_REQUIRE_EQUAL(edge.nodeIndex1, 0u);
      BOOST_REQUIRE_EQUAL(edge.nodeIndex2, 2u);
    }
    BOOST_REQUIRE(edger.next());

    edge = edger.getEdge();
    if (binIndex == 0) {
      BOOST_REQUIRE_EQUAL(edge.locusIndex, 0u);
      BOOST_REQUIRE_EQUAL(edge.nodeIndex1, 0u);
      BOOST_REQUIRE_EQUAL(edge.nodeIndex2, 3u);
    } else {
      BOOST_REQUIRE_EQUAL(edge.locusIndex, 1u);
      BOOST_REQUIRE_EQUAL(edge.nodeIndex1, 0u);
      BOOST_REQUIRE_EQUAL(edge.nodeIndex2, 1u);
    }

    BOOST_REQUIRE(!edger.next());
  }
}

BOOST_AUTO_TEST_CASE(test_EdgeRetrieverJumpManyBin2)
{
  SVLocus locus1;
  locusAddPair(locus1, 1, 10, 20, 2, 30, 40);
  SVLocus locus2;
  locusAddPair(locus2, 3, 10, 20, 4, 30, 40);
  SVLocus locus3;
  locusAddPair(locus3, 5, 10, 20, 6, 30, 40);
  SVLocus locus4;
  locusAddPair(locus4, 7, 10, 20, 8, 30, 40);
  SVLocus locus5;
  locusAddPair(locus5, 9, 10, 20, 10, 30, 40);
  SVLocus locus6;
  locusAddPair(locus6, 11, 10, 20, 12, 30, 40);

  SVLocusSetOptions sopt;
  sopt.minMergeEdgeObservations = 1;
  SVLocusSet set1(sopt);
  set1.merge(locus1);
  set1.merge(locus2);
  set1.merge(locus3);
  set1.merge(locus4);
  set1.merge(locus5);
  set1.merge(locus6);
  set1.checkState(true, true);

  static const unsigned binTotal(3);
  for (unsigned binIndex(0); binIndex < binTotal; ++binIndex) {
    EdgeRetrieverJumpBin edger(set1, 0, binTotal, binIndex);

    BOOST_REQUIRE(edger.next());

    EdgeInfo edge = edger.getEdge();
    BOOST_REQUIRE_EQUAL(edge.locusIndex, 0u + binIndex);
    BOOST_REQUIRE_EQUAL(edge.nodeIndex1, 0u);
    BOOST_REQUIRE_EQUAL(edge.nodeIndex2, 1u);

    BOOST_REQUIRE(edger.next());

    edge = edger.getEdge();
    BOOST_REQUIRE_EQUAL(edge.locusIndex, binTotal + binIndex);
    BOOST_REQUIRE_EQUAL(edge.nodeIndex1, 0u);
    BOOST_REQUIRE_EQUAL(edge.nodeIndex2, 1u);

    BOOST_REQUIRE(!edger.next());
  }
}

BOOST_AUTO_TEST_CASE(test_EdgeRetrieverJumpOddBin)
{
  SVLocus locus1;
  locusAddPair(locus1, 1, 10, 20, 2, 30, 40);
  SVLocus locus2;
  locusAddPair(locus2, 3, 10, 20, 4, 30, 40);
  SVLocus locus3;
  locusAddPair(locus3, 5, 10, 20, 6, 30, 40);
  SVLocus locus4;
  locusAddPair(locus4, 7, 10, 20, 8, 30, 40);
  SVLocus locus5;
  locusAddPair(locus5, 9, 10, 20, 10, 30, 40);
  SVLocus locus6;
  locusAddPair(locus6, 11, 10, 20, 12, 30, 40);
  SVLocus locus7;
  locusAddPair(locus7, 13, 10, 20, 14, 30, 40);

  SVLocusSetOptions sopt;
  sopt.minMergeEdgeObservations = 1;
  SVLocusSet set1(sopt);
  set1.merge(locus1);
  set1.merge(locus2);
  set1.merge(locus3);
  set1.merge(locus4);
  set1.merge(locus5);
  set1.merge(locus6);
  set1.merge(locus7);
  set1.checkState(true, true);

  unsigned              count(0);
  static const unsigned binTotal(3);
  for (unsigned binIndex(0); binIndex < binTotal; ++binIndex) {
    EdgeRetrieverJumpBin edger(set1, 0, binTotal, binIndex);

    while (edger.next()) {
      count++;
    }
  }

  BOOST_REQUIRE_EQUAL(count, 7u);
}

BOOST_AUTO_TEST_CASE(test_EdgeRetrieverJumpOddBinSelfEdge)
{
  SVLocus locus1;
  locusAddPair(locus1, 1, 10, 20, 1, 10, 20, true);
  SVLocus locus2;
  locusAddPair(locus2, 3, 10, 20, 3, 10, 20, true);
  SVLocus locus3;
  locusAddPair(locus3, 5, 10, 20, 5, 10, 20, true);
  SVLocus locus4;
  locusAddPair(locus4, 7, 10, 20, 8, 30, 40);
  SVLocus locus5;
  locusAddPair(locus5, 9, 10, 20, 9, 10, 20, true);
  SVLocus locus6;
  locusAddPair(locus6, 11, 10, 20, 12, 30, 40);
  SVLocus locus7;
  locusAddPair(locus7, 13, 10, 20, 14, 30, 40);

  locus1.mergeSelfOverlap();
  locus2.mergeSelfOverlap();
  locus3.mergeSelfOverlap();
  locus5.mergeSelfOverlap();

  SVLocusSetOptions sopt;
  sopt.minMergeEdgeObservations = 1;
  SVLocusSet set1(sopt);
  set1.merge(locus1);
  set1.merge(locus2);
  set1.merge(locus3);
  set1.merge(locus4);
  set1.merge(locus5);
  set1.merge(locus6);
  set1.merge(locus7);
  set1.checkState(true, true);

  std::set<unsigned> loci;
  unsigned           count(0);

  static const unsigned binTotal(3);
  for (unsigned binIndex(0); binIndex < binTotal; ++binIndex) {
    EdgeRetrieverJumpBin edger(set1, 0, binTotal, binIndex);

    while (edger.next()) {
      const EdgeInfo& edge(edger.getEdge());
      BOOST_REQUIRE_EQUAL(loci.count(edge.locusIndex), 0u);
      loci.insert(edge.locusIndex);
      count++;
    }
  }

  BOOST_REQUIRE_EQUAL(count, 7u);
}

BOOST_AUTO_TEST_SUITE_END()
