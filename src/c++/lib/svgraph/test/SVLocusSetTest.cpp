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

#include "svgraph/SVLocusSet.hpp"
#include "test/testSVLocusSetUtil.hpp"
#include "test/testSVLocusUtil.hpp"

#include "boost/timer/timer.hpp"

#include <fstream>  // For FindStringInFile

/// \brief Test the size and count of the properties of the SVLocusSet.
///
/// \param set1 The SV locus set to be tested
/// \param setSize The expected number of SVLocus objects in the set
/// \param setNonEmptySize The expected number of SVLocus objects that are not empty in the set
/// \param setNodeCount The expected total number of nodes in the SVLocus objects in the set
/// \param setEdgeCount The expected total number of edges in the SVLocus objects in the set
inline void TestSVLocusSetProperties(
    const SVLocusSet& set1, int setSize, int setNonEmptySize, int setNodeCount, int setEdgeCount)
{
  // Test the base SVLocusSet Contents
  BOOST_REQUIRE_EQUAL(set1.size(), setSize);
  BOOST_REQUIRE_EQUAL(set1.nonEmptySize(), setNonEmptySize);
  BOOST_REQUIRE_EQUAL(set1.totalNodeCount(), setNodeCount);
  BOOST_REQUIRE_EQUAL(set1.totalEdgeCount(), setEdgeCount);
}

/// \brief Test the distribution of the total edges and node observations for all nodes in the SVLocusSet.
///
/// \param set1 The SV locus set to be tested
/// \param edgeCountDistro The expected number of edges in each SVLocus object in the set
/// \param obsCountDistro The expected number of node observations in each SVLocus object in the set
static void TestSVLocusSetDistro(
    const SVLocusSet&            set1,
    const std::vector<unsigned>& edgeCountDistro,
    const std::vector<unsigned>& obsCountDistro)
{
  // review merged nodes
  static const unsigned maxCount(set1.totalEdgeCount());
  std::vector<unsigned> edgeCount(maxCount);
  std::vector<unsigned> obsCount(maxCount);

  set1.getNodeEdgeCountDistro(edgeCount);
  set1.getNodeObsCountDistro(obsCount);

  for (unsigned count = 0; count < maxCount; count++) {
    BOOST_REQUIRE_EQUAL(edgeCountDistro[count], edgeCount[count]);
  }

  for (unsigned count = 0; count < maxCount; count++) {
    BOOST_REQUIRE_EQUAL(obsCountDistro[count], obsCount[count]);
  }
}

/// \brief Test the distribution of out edges for a SV locus
///
/// \param set1 The SV locus set to be tested
/// \param locusIndex The locus index to be tested
/// \param expectedEdgeOutCount The expected number of out edges for each node of the tested locus
static void TestSVLocusSetOutEdges(
    const SVLocusSet& set1, unsigned locusIndex, const std::vector<unsigned>& expectedEdgeOutCount)
{
  assert(set1.getLocus(locusIndex).size() == expectedEdgeOutCount.size());

  for (unsigned count = 0; count < expectedEdgeOutCount.size(); count++) {
    BOOST_REQUIRE_EQUAL(set1.getLocus(locusIndex).getNode(count).outCount(), expectedEdgeOutCount[count]);
  }
}

/// \brief Check if a query string exists in a file name.
static bool FindStringInFile(const std::string& fileName, const std::string& searchString)
{
  std::ifstream iss(fileName);
  std::string   line;
  while (std::getline(iss, line)) {
    if (line.find(searchString) == std::string::npos) continue;
    return true;
  }
  return false;
}

BOOST_AUTO_TEST_SUITE(SVLocusSet_test_suite)

BOOST_AUTO_TEST_CASE(test_SVLocusSetProperties)
{
  // Create a SVLocusSet object and run it through the properties tests.
  SVLocusSetOptions sopt;
  sopt.minMergeEdgeObservations = 1;
  SVLocusSet set1(sopt);

  // Test getSource(..)
  BOOST_REQUIRE_EQUAL("UNKNOWN", set1.getSource());

  // Empty Set
  TestSVLocusSetProperties(set1, 0, 0, 0, 0);

  // Test merging empty locus
  SVLocus locus1;

  // Test empty(..)
  BOOST_REQUIRE(set1.empty());

  // Add a single node to the locus
  locus1.addNode(GenomeInterval(1, 10, 20));
  set1.merge(locus1);
  TestSVLocusSetProperties(set1, 1, 1, 1, 0);

  BOOST_REQUIRE(!set1.empty());

  // Add a second locus with 2 nodes linked by an edge.
  SVLocusSet          set2(sopt);
  SVLocus             locus2;
  const NodeIndexType nodePtr2 = locus2.addNode(GenomeInterval(2, 30, 40));
  const NodeIndexType nodePtr3 = locus2.addNode(GenomeInterval(3, 50, 60));
  locus2.linkNodes(nodePtr2, nodePtr3, 1, 0);
  set2.merge(locus2);
  TestSVLocusSetProperties(set2, 1, 1, 2, 2);

  SVLocusSet set3(sopt);
  set3.merge(locus1);
  set3.merge(locus2);

  TestSVLocusSetProperties(set3, 2, 2, 3, 2);

  BOOST_REQUIRE_EQUAL(set1.getAllSampleReadCounts().size(), 0);

  // Test setBuildTIme and addMergeTime. Can only test function does not fail.
  boost::timer::cpu_times tempTimes = {1, 2, 3};

  const CpuTimes cpuTimes(tempTimes);
  BOOST_REQUIRE_NO_THROW(set1.addBuildTime(cpuTimes));
  BOOST_REQUIRE_NO_THROW(set1.addMergeTime(cpuTimes));
}

BOOST_AUTO_TEST_CASE(test_SVLocusMergeLoci)
{
  BOOST_TEST_MESSAGE("SDS MANTA-699");

  // test merge of overlapping loci
  GenomeInterval testInterval1(1, 10, 20);
  GenomeInterval testInterval2(2, 30, 40);

  // construct a simple two-node locus
  SVLocus locus1;
  locusAddPair(locus1, 1, 10, 20, 2, 30, 40);

  SVLocus locus2;
  locusAddPair(locus2, 1, 10, 20, 2, 30, 40);

  SVLocusSetOptions sopt;
  sopt.minMergeEdgeObservations = 2;
  SVLocusSet        set1(sopt);
  const SVLocusSet& cset(set1);
  TestSVLocusSetProperties(cset, 0, 0, 0, 0);

  set1.merge(locus1);

  // Test the SVLocusSet with 1 node
  TestSVLocusSetProperties(cset, 1, 1, 2, 2);

  set1.merge(locus2);

  set1.checkState(true, true);

  // Test the SVLocusSet after merge
  TestSVLocusSetProperties(cset, 2, 1, 2, 2);

  bool testInterval1Found(false);
  bool testInterval2Found(false);
  for (const SVLocusNode& node : cset.getLocus(0)) {
    if (node.getInterval() == testInterval1) testInterval1Found = true;
    if (node.getInterval() == testInterval2) testInterval2Found = true;
  }

  BOOST_REQUIRE(testInterval1Found && testInterval2Found);
}

BOOST_AUTO_TEST_CASE(test_SVLocusSetMaxSearchCount)
{
  BOOST_TEST_MESSAGE("SDS MANTA-699");

  // Test that isUsable in merge fails when maxSearchCount is exceeded.
  GenomeInterval testInterval1(1, 10, 20);
  GenomeInterval testInterval2(2, 30, 40);

  // construct a simple two-node locus
  SVLocus locus1;
  locusAddPair(locus1, 1, 10, 20, 2, 30, 40);

  SVLocus locus2;
  locusAddPair(locus2, 1, 10, 20, 2, 30, 40);

  SVLocusSetOptions sopt;
  sopt.minMergeEdgeObservations = 2;
  sopt.maxSearchCount           = 1;
  SVLocusSet        set1(sopt);
  const SVLocusSet& cset(set1);
  TestSVLocusSetProperties(cset, 0, 0, 0, 0);

  set1.merge(locus1);
  set1.merge(locus2);

  set1.checkState(true, true);

  bool testInterval1Found(false);
  bool testInterval2Found(false);
  for (const SVLocusNode& node : cset.getLocus(0)) {
    if (node.getInterval() == testInterval1) testInterval1Found = true;
    if (node.getInterval() == testInterval2) testInterval2Found = true;
  }

  BOOST_REQUIRE(!(testInterval1Found || testInterval2Found));
}

BOOST_AUTO_TEST_CASE(test_SVLocusSetMaxSearchDensity)
{
  BOOST_TEST_MESSAGE("SDS MANTA-699");

  // Test that isUsable in merge fails when maxSearchDensity is exceeded.
  GenomeInterval testInterval1(1, 10, 20);
  GenomeInterval testInterval2(2, 30, 40);

  // construct a simple two-node locus
  SVLocus locus1;
  locusAddPair(locus1, 1, 10, 20, 2, 30, 40, false, 10000);

  SVLocus locus2;
  locusAddPair(locus2, 1, 10, 20, 2, 30, 40, false, 10000);

  SVLocusSetOptions sopt;
  sopt.minMergeEdgeObservations = 2;
  sopt.maxSearchDensity         = 0.01;
  SVLocusSet        set1(sopt);
  const SVLocusSet& cset(set1);
  TestSVLocusSetProperties(cset, 0, 0, 0, 0);

  set1.merge(locus1);
  set1.merge(locus2);

  set1.checkState(true, true);

  bool testInterval1Found(false);
  bool testInterval2Found(false);
  for (const SVLocusNode& node : cset.getLocus(0)) {
    if (node.getInterval() == testInterval1) testInterval1Found = true;
    if (node.getInterval() == testInterval2) testInterval2Found = true;
  }

  BOOST_REQUIRE(!(testInterval1Found || testInterval2Found));
}

BOOST_AUTO_TEST_CASE(test_SVLocusMultiOverlapMerge)
{
  BOOST_TEST_MESSAGE("SDS MANTA-699");

  // test merge of overlapping loci, reproduces production failure
  SVLocus locus1;
  locusAddPair(locus1, 1, 10, 20, 12, 30, 40);

  SVLocus locus2;
  locusAddPair(locus2, 2, 10, 20, 12, 50, 60);

  SVLocus locus3;
  locusAddPair(locus3, 3, 10, 20, 12, 35, 55);

  SVLocusSetOptions sopt;
  sopt.minMergeEdgeObservations = 1;
  SVLocusSet set1(sopt);
  set1.merge(locus1);
  set1.merge(locus2);
  set1.merge(locus3);
  set1.checkState(true, true);
  const SVLocusSet& cset1(set1);

  GenomeInterval testInterval(12, 30, 60);

  BOOST_REQUIRE_EQUAL(cset1.nonEmptySize(), 1u);
  BOOST_REQUIRE_EQUAL(cset1.getLocus(0).size(), 4u);

  bool isFound(false);
  for (const SVLocusNode& node : cset1.getLocus(0)) {
    if (node.getInterval() == testInterval) isFound = true;
  }
  BOOST_REQUIRE(isFound);

  // review merged nodes
  const std::vector<unsigned> expectedEdgeCount = {0, 3, 0, 1, 0, 0};
  const std::vector<unsigned> expectedObsCount  = {1, 3, 0, 0, 0, 0};

  TestSVLocusSetDistro(set1, expectedEdgeCount, expectedObsCount);

  const std::vector<unsigned> expectedOutEdgeCount = {1, 0, 1, 1};
  TestSVLocusSetOutEdges(set1, 0, expectedOutEdgeCount);
}

BOOST_AUTO_TEST_CASE(test_SVLocusMultiOverlapMerge2)
{
  // test merge of overlapping loci, reproduces production failure

  SVLocus locus1;
  {
    NodeIndexType nodePtr1 = locus1.addNode(GenomeInterval(1, 10, 20));
    NodeIndexType nodePtr2 = locus1.addNode(GenomeInterval(1, 30, 40));
    NodeIndexType nodePtr3 = locus1.addNode(GenomeInterval(1, 50, 60));
    locus1.linkNodes(nodePtr1, nodePtr2);
    locus1.linkNodes(nodePtr1, nodePtr3);
  }

  SVLocus locus2;
  locusAddPair(locus2, 1, 10, 60, 2, 10, 60);

  SVLocusSetOptions sopt;
  sopt.minMergeEdgeObservations = 1;
  SVLocusSet set1(sopt);
  set1.merge(locus1);
  set1.merge(locus2);
  set1.checkState(true, true);
  const SVLocusSet& cset1(set1);

  GenomeInterval testInterval(1, 10, 60);

  BOOST_REQUIRE_EQUAL(cset1.nonEmptySize(), 1u);
  BOOST_REQUIRE_EQUAL(cset1.getLocus(0).size(), 2u);

  bool isFound(false);
  for (const SVLocusNode& node : cset1.getLocus(0)) {
    if (node.getInterval() == testInterval) isFound = true;
  }
  BOOST_REQUIRE(isFound);

  // review merged nodes
  BOOST_REQUIRE_EQUAL(cset1.totalEdgeCount(), 3);

  // nodePtr1-3 merge into a single node 10-60, matching the testInverval
  // Because first node of testInverval is id 1, 3 outgoing edges, and 0 for other.
  const std::vector<unsigned> expectedEdgeCount = {0, 1, 1, 0, 0, 0};
  const std::vector<unsigned> expectedObsCount  = {1, 0, 0, 1, 0, 0};

  TestSVLocusSetDistro(set1, expectedEdgeCount, expectedObsCount);
}

BOOST_AUTO_TEST_CASE(test_SVLocusMultiOverlapMerge3)
{
  // test merge of overlapping loci, reproduces production failure

  SVLocus locus1;
  locusAddPair(locus1, 1, 10, 20, 3, 10, 20);

  SVLocus locus2;
  locusAddPair(locus2, 1, 30, 40, 4, 10, 20);

  SVLocus locus3;
  locusAddPair(locus3, 2, 30, 40, 5, 10, 20);

  SVLocus locus4;
  locusAddPair(locus4, 1, 15, 35, 6, 10, 20);

  SVLocus locus5;
  locusAddPair(locus5, 2, 15, 35, 7, 10, 20);

  SVLocusSetOptions sopt;
  sopt.minMergeEdgeObservations = 1;
  SVLocusSet set1(sopt);
  set1.merge(locus1);
  set1.merge(locus2);
  set1.merge(locus3);
  set1.merge(locus4);
  set1.merge(locus5);
  set1.checkState(true, true);
  const SVLocusSet& cset1(set1);

  GenomeInterval testInterval(1, 10, 40);

  BOOST_REQUIRE_EQUAL(cset1.nonEmptySize(), 2u);
  BOOST_REQUIRE_EQUAL(cset1.getLocus(0).size(), 4u);

  bool isFound(false);
  for (const SVLocusNode& node : cset1.getLocus(0)) {
    if (node.getInterval() == testInterval) isFound = true;
  }
  BOOST_REQUIRE(isFound);

  // review merged nodes
  BOOST_REQUIRE_EQUAL(cset1.totalEdgeCount(), 10);

  // nodePtr1-3 merge into a single node 10-60, matching the testInverval
  // Because first node of testInverval is id 1, 3 outgoing edges, and 0 for other.
  const std::vector<unsigned> expectedEdgeCount = {0, 5, 1, 1, 0, 0};
  const std::vector<unsigned> expectedObsCount  = {5, 0, 1, 1, 0, 0};

  TestSVLocusSetDistro(cset1, expectedEdgeCount, expectedObsCount);
}

BOOST_AUTO_TEST_CASE(test_SVLocusMultiOverlapMerge4)
{
  // test merge of overlapping loci, reproduces production failure

  SVLocus locus1;
  locusAddPair(locus1, 1, 10, 60, 2, 20, 30);

  SVLocus locus2;
  locusAddPair(locus2, 1, 40, 50, 1, 20, 30);

  SVLocusSetOptions sopt;
  sopt.minMergeEdgeObservations = 1;
  SVLocusSet set1(sopt);
  set1.merge(locus1);
  set1.merge(locus2);

  const SVLocusSet& cset1(set1);
  cset1.checkState(true, true);

  GenomeInterval testInterval(1, 10, 60);

  BOOST_REQUIRE_EQUAL(cset1.nonEmptySize(), 1u);
  BOOST_REQUIRE_EQUAL(cset1.getLocus(0).size(), 2u);

  bool isFound(false);
  for (const SVLocusNode& node : cset1.getLocus(0)) {
    if (node.getInterval() == testInterval) isFound = true;
  }
  BOOST_REQUIRE(isFound);

  // review merged nodes
  BOOST_REQUIRE_EQUAL(cset1.totalEdgeCount(), 3);

  const std::vector<unsigned> expectedEdgeCount = {0, 1, 1, 0, 0, 0};
  const std::vector<unsigned> expectedObsCount  = {1, 0, 1, 0, 0, 0};

  TestSVLocusSetDistro(cset1, expectedEdgeCount, expectedObsCount);
}

BOOST_AUTO_TEST_CASE(test_SVLocusNoiseMerge)
{
  BOOST_TEST_MESSAGE("SDS MANTA-699");
  SVLocus locus1;
  locusAddPair(locus1, 1, 10, 60, 2, 20, 30);

  SVLocus locus2;
  locusAddPair(locus2, 1, 10, 60, 2, 20, 30);

  SVLocus locus3;
  locusAddPair(locus3, 1, 10, 60, 3, 20, 30);

  // Test with min observations set to 1
  SVLocusSetOptions sopt1;
  sopt1.minMergeEdgeObservations = 1;
  SVLocusSet set1(sopt1);
  set1.merge(locus1);
  set1.merge(locus2);
  set1.merge(locus3);
  const SVLocusSet& cset1(set1);

  cset1.checkState(true, true);
  BOOST_REQUIRE_EQUAL(cset1.nonEmptySize(), 1u);
  BOOST_REQUIRE_EQUAL(cset1.getLocus(0).size(), 3u);

  // review merged nodes
  BOOST_REQUIRE_EQUAL(cset1.totalEdgeCount(), 4);
  const std::vector<unsigned> expectedEdgeCount1 = {0, 2, 1, 0, 0, 0};
  const std::vector<unsigned> expectedObsCount1  = {2, 0, 0, 1, 0, 0};
  TestSVLocusSetDistro(cset1, expectedEdgeCount1, expectedObsCount1);

  // Test with min observations set to 2
  SVLocusSetOptions sopt2;
  sopt2.minMergeEdgeObservations = 2;
  SVLocusSet set2(sopt2);
  set2.merge(locus1);
  set2.merge(locus2);
  set2.merge(locus3);
  const SVLocusSet& cset2(set2);

  cset2.checkState(true, true);
  BOOST_REQUIRE_EQUAL(cset2.nonEmptySize(), 2u);
  BOOST_REQUIRE_EQUAL(cset2.getLocus(0).size(), 2u);

  // review merged nodes
  BOOST_REQUIRE_EQUAL(cset2.totalEdgeCount(), 4);
  const std::vector<unsigned> expectedEdgeCount2 = {0, 4, 0, 0, 0, 0};
  const std::vector<unsigned> expectedObsCount2  = {2, 1, 1, 0, 0, 0};
  TestSVLocusSetDistro(cset2, expectedEdgeCount2, expectedObsCount2);

  // Test with no min observations
  SVLocusSetOptions sopt3;
  sopt3.minMergeEdgeObservations = 3;
  SVLocusSet set3(sopt3);
  set3.merge(locus1);
  set3.merge(locus2);
  set3.merge(locus3);
  const SVLocusSet& cset3(set3);

  cset3.checkState(true, true);
  BOOST_REQUIRE_EQUAL(cset3.nonEmptySize(), 3u);
  BOOST_REQUIRE_EQUAL(cset3.getLocus(0).size(), 2u);

  // review merged nodes
  BOOST_REQUIRE_EQUAL(cset3.totalEdgeCount(), 6);
  const std::vector<unsigned> expectedEdgeCount3 = {0, 6, 0, 0, 0, 0};
  const std::vector<unsigned> expectedObsCount3  = {3, 3, 0, 0, 0, 0};
  TestSVLocusSetDistro(cset3, expectedEdgeCount3, expectedObsCount3);
}

BOOST_AUTO_TEST_CASE(test_SVLocusNoiseClean)
{
  SVLocus locus1;
  locusAddPair(locus1, 1, 10, 60, 2, 20, 30);

  SVLocus locus2;
  locusAddPair(locus2, 1, 10, 60, 2, 20, 30);

  SVLocus locus3;
  locusAddPair(locus3, 1, 10, 60, 3, 20, 30);

  {
    // clean out any nodes that are below the merge threshold.
    SVLocusSetOptions sopt;
    sopt.minMergeEdgeObservations = 2;
    SVLocusSet set1(sopt);
    set1.merge(locus1);
    set1.merge(locus2);
    set1.merge(locus3);
    const SVLocusSet& cset1(set1);

    cset1.checkState(true, true);
    TestSVLocusSetProperties(cset1, 2, 2, 4, 4);

    set1.clean();
    TestSVLocusSetProperties(cset1, 2, 1, 2, 2);
  }

  {
    // clean out any nodes in each GenomeInterval that are below the merge threshold.
    SVLocusSetOptions sopt;
    sopt.minMergeEdgeObservations = 2;
    SVLocusSet set1(sopt);
    set1.merge(locus1);
    set1.merge(locus2);
    set1.merge(locus3);
    const SVLocusSet& cset1(set1);

    TestSVLocusSetProperties(cset1, 2, 2, 4, 4);
    BOOST_REQUIRE_EQUAL(cset1.getLocus(1).size(), 2u);

    set1.cleanRegion(GenomeInterval(3, 0, 70));
    TestSVLocusSetProperties(cset1, 3, 2, 4, 4);
    BOOST_REQUIRE_EQUAL(cset1.getLocus(1).size(), 2u);

    set1.cleanRegion(GenomeInterval(1, 0, 70));
    TestSVLocusSetProperties(cset1, 3, 1, 2, 2);
    BOOST_REQUIRE_EQUAL(cset1.getLocus(0).size(), 2u);
  }
}

BOOST_AUTO_TEST_CASE(test_SVLocusNoiseCleanOrder)
{
  SVLocus locus1;
  locusAddPair(locus1, 1, 10, 60, 2, 20, 30);

  SVLocus locus2;
  locusAddPair(locus2, 1, 10, 60, 2, 20, 30);

  SVLocus locus3;
  locusAddPair(locus3, 1, 10, 60, 3, 20, 30);

  SVLocus locus4;
  locusAddPair(locus4, 1, 10, 60, 4, 20, 30);

  SVLocus locus5;
  locusAddPair(locus5, 1, 10, 60, 4, 20, 30);

  SVLocus locus6;
  locusAddPair(locus6, 1, 10, 60, 5, 20, 30);

  {
    SVLocusSetOptions sopt;
    sopt.minMergeEdgeObservations = 2;
    SVLocusSet set1(sopt);
    set1.merge(locus1);
    set1.merge(locus2);
    set1.merge(locus3);
    set1.merge(locus4);
    set1.merge(locus5);
    set1.merge(locus6);
    const SVLocusSet& cset1(set1);

    TestSVLocusSetProperties(cset1, 4, 3, 7, 8);
    BOOST_REQUIRE_EQUAL(cset1.getLocus(1).size(), 2u);
    BOOST_REQUIRE_EQUAL(cset1.getLocus(0).getNode(0).getInterval().tid, 1);
    BOOST_REQUIRE_EQUAL(cset1.getLocus(0).getNode(1).getInterval().tid, 2);
    BOOST_REQUIRE_EQUAL(cset1.getLocus(0).getNode(2).getInterval().tid, 4);
    BOOST_REQUIRE_EQUAL(cset1.getLocus(1).getNode(1).getInterval().tid, 3);
    BOOST_REQUIRE_EQUAL(cset1.getLocus(2).getNode(1).getInterval().tid, 5);

    set1.cleanRegion(GenomeInterval(1, 0, 70));

    TestSVLocusSetProperties(cset1, 4, 1, 3, 4);
    BOOST_REQUIRE_EQUAL(cset1.getLocus(0).size(), 3u);
    BOOST_REQUIRE_EQUAL(cset1.getLocus(0).getNode(0).getInterval().tid, 1);
    BOOST_REQUIRE_EQUAL(cset1.getLocus(0).getNode(1).getInterval().tid, 2);
    BOOST_REQUIRE_EQUAL(cset1.getLocus(0).getNode(2).getInterval().tid, 4);
  }
}

BOOST_AUTO_TEST_CASE(test_SVLocusNoiseCleanRemote)
{
  /// locus1 has 2 nodes of the same id linked. cleanRegion will remove it.
  SVLocus locus1;
  locusAddPair(locus1, 1, 100, 110, 1, 10, 20);

  {
    SVLocusSetOptions sopt;
    sopt.minMergeEdgeObservations = 2;
    SVLocusSet set1(sopt);
    set1.merge(locus1);
    const SVLocusSet& cset1(set1);

    TestSVLocusSetProperties(cset1, 1, 1, 2, 2);
    BOOST_REQUIRE_EQUAL(cset1.getLocus(0).size(), 2u);

    set1.cleanRegion(GenomeInterval(1, 0, 120));
    TestSVLocusSetProperties(cset1, 2, 0, 0, 0);
    BOOST_REQUIRE_EQUAL(cset1.getLocus(0).size(), 0u);
  }
}

BOOST_AUTO_TEST_CASE(test_SVLocusEvidenceRange)
{
  BOOST_TEST_MESSAGE("SDS MANTA-699");

  /// locus 1 and 2 will be merged and the evidence range for the merged locus
  /// will enlarge to the range of both of the known_pos_range2.
  SVLocus locus1;
  {
    NodeIndexType node1 = locus1.addNode(GenomeInterval(1, 100, 110));
    NodeIndexType node2 = locus1.addNode(GenomeInterval(2, 100, 110));
    locus1.linkNodes(node1, node2);
    locus1.setNodeEvidence(node1, known_pos_range2(50, 60));
  }

  SVLocus locus2;
  {
    NodeIndexType node1 = locus2.addNode(GenomeInterval(1, 100, 110));
    NodeIndexType node2 = locus2.addNode(GenomeInterval(2, 100, 110));
    locus2.linkNodes(node1, node2);
    locus2.setNodeEvidence(node1, known_pos_range2(30, 40));
  }

  {
    SVLocusSetOptions sopt;
    sopt.minMergeEdgeObservations = 2;
    SVLocusSet set1(sopt);
    set1.merge(locus1);
    set1.merge(locus2);
    const SVLocusSet& cset1(set1);

    TestSVLocusSetProperties(cset1, 2, 1, 2, 2);
    BOOST_REQUIRE_EQUAL(cset1.getLocus(0).getNode(0).getEvidenceRange(), known_pos_range2(30, 60));
  }
}

BOOST_AUTO_TEST_CASE(test_SVLocusNoiseOverlap)
{
  BOOST_TEST_MESSAGE("SDS MANTA-700");

  // adapted from a production failure case:
  /// 4 loci, locus1,2 overlap perfectly, locus3 partially overlaps
  /// locus1,2 and locus 4 overlaps locus 3.
  /// locus3,4 2nd nodes are the same.
  /// All pairs merge into a single locus with locus0 range 10-75.
  SVLocus locus1;
  locusAddPair(locus1, 1, 10, 60, 2, 20, 30);
  SVLocus locus2;
  locusAddPair(locus2, 1, 10, 60, 2, 20, 30);
  SVLocus locus3;
  locusAddPair(locus3, 1, 59, 70, 3, 20, 30);
  SVLocus locus4;
  locusAddPair(locus4, 1, 65, 75, 3, 20, 30);

  {
    SVLocusSetOptions sopt;
    sopt.minMergeEdgeObservations = 2;
    SVLocusSet set1(sopt);
    set1.merge(locus1);
    set1.merge(locus2);
    set1.merge(locus3);
    set1.merge(locus4);
    const SVLocusSet& cset1(set1);

    set1.finalize();
    cset1.checkState(true, true);

    TestSVLocusSetProperties(cset1, 3, 1, 3, 4);
    BOOST_REQUIRE_EQUAL(cset1.getLocus(0).size(), 3u);
    BOOST_REQUIRE_EQUAL(cset1.getLocus(0).getNode(0).getInterval(), GenomeInterval(1, 10, 75));
  }
}

BOOST_AUTO_TEST_CASE(test_SVLocusSingleSelfEdge)
{
  BOOST_TEST_MESSAGE("SDS MANTA-699");

  SVLocus locus1;
  locusAddPair(locus1, 1, 10, 60, 1, 20, 70);
  locus1.mergeSelfOverlap();

  {
    SVLocusSetOptions sopt;
    sopt.minMergeEdgeObservations = 1;
    SVLocusSet set1(sopt);
    set1.merge(locus1);
    const SVLocusSet& cset1(set1);

    cset1.checkState(true, true);
    TestSVLocusSetProperties(cset1, 1, 1, 1, 1);
    BOOST_REQUIRE_EQUAL(cset1.getLocus(0).size(), 1u);
    BOOST_REQUIRE_EQUAL(cset1.selfEdgeCount(), 1);

    set1.clean();

    TestSVLocusSetProperties(cset1, 1, 1, 1, 1);
    BOOST_REQUIRE_EQUAL(cset1.getLocus(0).size(), 1u);
  }
}

BOOST_AUTO_TEST_CASE(test_SVLocusDoubleSelfEdge)
{
  BOOST_TEST_MESSAGE("SDS MANTA-699");

  /// 2 loci have self edges, locus2 is set to both local nodes.
  /// After merging the locus has a single node pointing to itself twice.
  SVLocus locus1;
  locusAddPair(locus1, 1, 10, 60, 1, 20, 70);
  SVLocus locus2;
  locusAddPair(locus2, 1, 20, 70, 1, 10, 60, true);

  locus1.mergeSelfOverlap();
  locus2.mergeSelfOverlap();

  {
    SVLocusSetOptions sopt;
    sopt.minMergeEdgeObservations = 1;
    SVLocusSet set1(sopt);
    set1.merge(locus1);
    set1.merge(locus2);
    const SVLocusSet& cset1(set1);

    cset1.checkState(true, true);
    TestSVLocusSetProperties(cset1, 2, 1, 1, 1);
    BOOST_REQUIRE_EQUAL(cset1.getLocus(0).size(), 1u);
    BOOST_REQUIRE_EQUAL(cset1.getLocus(0).getNode(0).size(), 1u);
    BOOST_REQUIRE_EQUAL(cset1.getLocus(0).getNode(0).getInterval(), GenomeInterval(1, 10, 70));
    BOOST_REQUIRE_EQUAL(cset1.getLocus(0).getNode(0).outCount(), 2u);
  }
}

BOOST_AUTO_TEST_CASE(test_SVLocusDoubleSelfEdge2)
{
  // Test merging of normal locus with double self edge locus.
  //
  // Locus1 is a default pair, Locus2 has a self pair.
  // Locus1 and Locus 2 merge. node1 will point at self and nodeId 2.
  SVLocus locus1;
  locusAddPair(locus1, 1, 10, 60, 2, 20, 70);
  SVLocus locus2;
  locusAddPair(locus2, 1, 10, 60, 1, 10, 60, true);

  locus1.mergeSelfOverlap();
  locus2.mergeSelfOverlap();

  {
    SVLocusSetOptions sopt;
    sopt.minMergeEdgeObservations = 1;
    SVLocusSet set1(sopt);
    set1.merge(locus1);
    set1.merge(locus2);
    const SVLocusSet& cset1(set1);

    cset1.checkState(true, true);
    TestSVLocusSetProperties(cset1, 2, 1, 2, 3);
    BOOST_REQUIRE_EQUAL(cset1.nonEmptySize(), 1u);
    BOOST_REQUIRE_EQUAL(cset1.getLocus(0).getNode(0).size(), 2u);
    BOOST_REQUIRE_EQUAL(cset1.getLocus(0).getNode(0).outCount(), 2u);
    BOOST_REQUIRE_EQUAL(cset1.getLocus(0).getNode(0).getInterval(), GenomeInterval(1, 10, 60));
    BOOST_REQUIRE_EQUAL(cset1.getLocus(0).getNode(1).getInterval(), GenomeInterval(2, 20, 70));
  }
}

BOOST_AUTO_TEST_CASE(test_SVLocusNodeOverlapEdge)
{
  BOOST_TEST_MESSAGE("SDS MANTA-699");

  // test merge of two edges: one edge node encompasses both nodes of the second edge:
  //
  // Locus2 nodes merge with locus1-node1.
  // Nothing merges because they do not reach the minMergedEdgeObservations.

  SVLocus locus1;
  locusAddPair(locus1, 1, 10, 60, 2, 20, 70);
  SVLocus locus2;
  locusAddPair(locus2, 1, 10, 20, 1, 30, 40);

  locus1.mergeSelfOverlap();
  locus2.mergeSelfOverlap();

  {
    SVLocusSetOptions sopt;
    sopt.minMergeEdgeObservations = 2;
    SVLocusSet set1(sopt);
    set1.merge(locus1);
    set1.merge(locus2);
    const SVLocusSet& cset1(set1);

    BOOST_REQUIRE_EQUAL(cset1.nonEmptySize(), 2u);
    set1.finalize();
    cset1.checkState(true, true);
    TestSVLocusSetProperties(cset1, 2, 0, 0, 0);
  }

  {
    // reverse the order of locus addition to be sure:
    SVLocusSetOptions sopt;
    sopt.minMergeEdgeObservations = 2;
    SVLocusSet set1(sopt);
    set1.merge(locus2);
    set1.merge(locus1);
    const SVLocusSet& cset1(set1);

    BOOST_REQUIRE_EQUAL(cset1.nonEmptySize(), 2u);

    set1.finalize();
    cset1.checkState(true, true);
    TestSVLocusSetProperties(cset1, 2, 0, 0, 0);
  }
}

BOOST_AUTO_TEST_CASE(test_SVLocusNodeOverlapSelfEdge)
{
  // test merge of two edges: one edge node overlaps a self-edge
  //
  // Locus1 is a standard locus, Locus2 is a self edge locus.
  // Nothing merges because they do not reach the minMergedEdgeObservations.
  //
  // Last, add locus3 self edge idential to locus 2, they merge, saved after finalize.

  SVLocus locus1;
  locusAddPair(locus1, 1, 10, 60, 2, 20, 70);
  SVLocus locus2;
  locusAddPair(locus2, 1, 10, 20, 1, 10, 20, true);

  locus1.mergeSelfOverlap();
  locus2.mergeSelfOverlap();

  {
    SVLocusSetOptions sopt;
    sopt.minMergeEdgeObservations = 2;
    SVLocusSet set1(sopt);
    set1.merge(locus1);
    set1.merge(locus2);
    const SVLocusSet& cset1(set1);

    BOOST_REQUIRE_EQUAL(cset1.nonEmptySize(), 2u);

    set1.finalize();
    cset1.checkState(true, true);
    TestSVLocusSetProperties(cset1, 2, 0, 0, 0);
  }

  {
    // reverse the order of locus addition to be sure:
    SVLocusSetOptions sopt;
    sopt.minMergeEdgeObservations = 2;
    SVLocusSet set1(sopt);
    set1.merge(locus2);
    set1.merge(locus1);
    const SVLocusSet& cset1(set1);

    BOOST_REQUIRE_EQUAL(cset1.nonEmptySize(), 2u);

    set1.finalize();
    cset1.checkState(true, true);
    BOOST_REQUIRE_EQUAL(cset1.nonEmptySize(), 0u);
    TestSVLocusSetProperties(cset1, 2, 0, 0, 0);
  }

  SVLocus locus3;
  locusAddPair(locus3, 1, 10, 20, 1, 10, 20, true);

  locus3.mergeSelfOverlap();

  {
    SVLocusSetOptions sopt;
    sopt.minMergeEdgeObservations = 2;
    SVLocusSet set1(sopt);
    set1.merge(locus1);
    set1.merge(locus2);
    set1.merge(locus3);
    const SVLocusSet& cset1(set1);

    BOOST_REQUIRE_EQUAL(cset1.nonEmptySize(), 2u);

    set1.finalize();
    cset1.checkState(true, true);
    BOOST_REQUIRE_EQUAL(cset1.nonEmptySize(), 1u);
    TestSVLocusSetProperties(cset1, 3, 1, 1, 1);
  }
}

BOOST_AUTO_TEST_CASE(test_SVLocusMergeToSelfEdge)
{
  BOOST_TEST_MESSAGE("SDS MANTA-699");

  // test merge or edges which are not self-edges themselves, but merge to
  // a state where they should be a self edge.
  //
  // Locus 1 and 2 merge into a single locus with interval 10-40 with 2 self edges.

  SVLocus locus1;
  locusAddPair(locus1, 1, 10, 20, 1, 25, 40);
  SVLocus locus2;
  locusAddPair(locus2, 1, 15, 30, 1, 35, 40);

  locus1.mergeSelfOverlap();
  locus2.mergeSelfOverlap();

  {
    SVLocusSetOptions sopt;
    sopt.minMergeEdgeObservations = 2;
    SVLocusSet set1(sopt);
    set1.merge(locus1);
    set1.merge(locus2);
    const SVLocusSet& cset1(set1);

    set1.finalize();
    cset1.checkState(true, true);
    BOOST_REQUIRE_EQUAL(cset1.nonEmptySize(), 1u);
    TestSVLocusSetProperties(cset1, 2, 1, 1, 1);
  }
}

BOOST_AUTO_TEST_CASE(test_SVLocusMergeToSelfEdge2)
{
  // Situation: Large signal region spans connected region pair.
  //
  // Criteria: Large region gains self-edge only if connected regions are signal

  SVLocus locus1;
  locusAddPair(locus1, 1, 10, 20, 1, 30, 40);
  SVLocus locus2;
  locusAddPair(locus2, 1, 10, 20, 1, 30, 40);
  SVLocus locus3;
  locusAddPair(locus3, 1, 10, 40, 2, 10, 20);
  SVLocus locus4;
  locusAddPair(locus4, 1, 10, 40, 2, 10, 20);

  {
    // test non-signal spanned pair (skipped locus2)
    // locus1 does not get merged with other loci due to low signal.
    SVLocusSetOptions sopt;
    sopt.minMergeEdgeObservations = 2;
    SVLocusSet set1(sopt);
    set1.merge(locus1);
    set1.merge(locus3);
    set1.merge(locus4);
    const SVLocusSet& cset1(set1);

    set1.finalize();
    cset1.checkState(true, true);
    TestSVLocusSetProperties(cset1, 3, 1, 2, 2);
    BOOST_REQUIRE_EQUAL(cset1.nonEmptySize(), 1u);
    BOOST_REQUIRE_EQUAL(cset1.getLocus(1).size(), 2u);
    BOOST_REQUIRE_EQUAL(cset1.getLocus(1).getNode(0).outCount(), 2);
    BOOST_REQUIRE_EQUAL(cset1.getLocus(1).getNode(0).size(), 1u);
  }

  {
    // test signal spanned pair
    // locus1 gets merged with other loci due to strong enough signal.

    SVLocusSetOptions sopt;
    sopt.minMergeEdgeObservations = 2;
    SVLocusSet set1(sopt);
    set1.merge(locus1);
    set1.merge(locus2);
    set1.merge(locus3);
    set1.merge(locus4);
    const SVLocusSet& cset1(set1);

    set1.finalize();
    cset1.checkState(true, true);
    TestSVLocusSetProperties(cset1, 3, 1, 2, 3);
    BOOST_REQUIRE_EQUAL(cset1.nonEmptySize(), 1u);
    BOOST_REQUIRE_EQUAL(cset1.getLocus(0).size(), 2u);
    BOOST_REQUIRE_EQUAL(cset1.getLocus(0).getNode(0).outCount(), 4);
    BOOST_REQUIRE_EQUAL(cset1.getLocus(0).getNode(0).size(), 2u);
  }
}

BOOST_AUTO_TEST_CASE(test_SVLocusMergeToSelfEdge3)
{
  // Situation: Large signal region with self edge spans connected region.
  //
  // Criteria: self edge count increases by one

  SVLocus locus1;
  locusAddPair(locus1, 1, 10, 40, 1, 10, 40, true);
  SVLocus locus2;
  locusAddPair(locus2, 1, 10, 40, 1, 10, 40, true);
  SVLocus locus3;
  locusAddPair(locus3, 1, 10, 20, 1, 30, 40);

  locus1.mergeSelfOverlap();
  locus2.mergeSelfOverlap();

  {
    SVLocusSetOptions sopt;
    sopt.minMergeEdgeObservations = 2;
    SVLocusSet set1(sopt);
    set1.merge(locus1);
    set1.merge(locus2);
    set1.merge(locus3);
    const SVLocusSet& cset1(set1);

    set1.finalize();
    cset1.checkState(true, true);
    TestSVLocusSetProperties(cset1, 2, 1, 1, 1);
    BOOST_REQUIRE_EQUAL(cset1.getLocus(0).size(), 1u);
    BOOST_REQUIRE_EQUAL(cset1.getLocus(0).getNode(0).size(), 1u);
    BOOST_REQUIRE_EQUAL(cset1.getLocus(0).getNode(0).outCount(), 3u);
  }

  {
    // run again with locus3 first
    SVLocusSetOptions sopt;
    sopt.minMergeEdgeObservations = 2;
    SVLocusSet set1(sopt);
    set1.merge(locus3);
    set1.merge(locus1);
    set1.merge(locus2);
    const SVLocusSet& cset1(set1);

    set1.finalize();
    cset1.checkState(true, true);
    TestSVLocusSetProperties(cset1, 2, 1, 1, 1);
    BOOST_REQUIRE_EQUAL(cset1.getLocus(0).size(), 1u);
    BOOST_REQUIRE_EQUAL(cset1.getLocus(0).getNode(0).size(), 1u);
    BOOST_REQUIRE_EQUAL(cset1.getLocus(0).getNode(0).outCount(), 3u);
  }
}

#if 0
BOOST_AUTO_TEST_CASE( test_SVLocusMergeToSelfEdge4 )
{
    // Situation: Large signal region spans 2 connected region pairs.
    //
    // Criteria: Large region gains self-edge only if the connected region pairs total to a self-edge signal


    SVLocus locus1;
    locusAddPair(locus1,1,10,20,1,30,40);
    SVLocus locus2;
    locusAddPair(locus2,1,10,20,1,30,40);
    SVLocus locus3;
    locusAddPair(locus3,1,10,40,2,10,20);
    SVLocus locus4;
    locusAddPair(locus4,1,10,40,2,10,20);

    {
        // test non-signal spanned pair
        SVLocusSet set1(2);
        set1.merge(locus1);
        set1.merge(locus3);
        set1.merge(locus4);
        const SVLocusSet& cset1(set1);

        set1.finalize();
        cset1.checkState(true,true);
        BOOST_REQUIRE_EQUAL(cset1.nonEmptySize(),1u);
        BOOST_REQUIRE_EQUAL(cset1.getLocus(1).size(),2u);
        BOOST_REQUIRE_EQUAL(cset1.getLocus(1).getNode(0).size(),1u);
    }

    {
        // test signal spanned pair
        SVLocusSet set1(2);
        set1.merge(locus1);
        set1.merge(locus2);
        set1.merge(locus3);
        set1.merge(locus4);
        const SVLocusSet& cset1(set1);

        set1.finalize();
        cset1.checkState(true,true);
        BOOST_REQUIRE_EQUAL(cset1.nonEmptySize(),1u);
        BOOST_REQUIRE_EQUAL(cset1.getLocus(0).size(),2u);
        BOOST_REQUIRE_EQUAL(cset1.getLocus(0).getNode(0).size(),2u);
    }
}
#endif

BOOST_AUTO_TEST_CASE(test_SVLocusSmallRegionClean)
{
  //
  // challenge the region cleaner with regions which (1) span both nodes of a pair (2) span the first node (3)
  // span the second node
  //
  // cleanRegion(..) will remove loci based on overlapping the node with the out edge.
  SVLocus locus1;
  locusAddPair(locus1, 1, 10, 20, 1, 30, 40);

  {
    SVLocusSetOptions sopt;
    sopt.minMergeEdgeObservations = 2;
    SVLocusSet set1(sopt);
    set1.merge(locus1);
    const SVLocusSet& cset1(set1);

    set1.cleanRegion(GenomeInterval(1, 0, 70));

    TestSVLocusSetProperties(cset1, 2, 0, 0, 0);
  }

  {
    SVLocusSetOptions sopt;
    sopt.minMergeEdgeObservations = 2;
    SVLocusSet set1(sopt);
    set1.merge(locus1);
    const SVLocusSet& cset1(set1);

    set1.cleanRegion(GenomeInterval(1, 25, 70));

    TestSVLocusSetProperties(cset1, 2, 1, 2, 2);
  }

  {
    SVLocusSetOptions sopt;
    sopt.minMergeEdgeObservations = 2;
    SVLocusSet set1(sopt);
    set1.merge(locus1);
    const SVLocusSet& cset1(set1);

    set1.cleanRegion(GenomeInterval(1, 5, 25));

    TestSVLocusSetProperties(cset1, 2, 0, 0, 0);
  }

  {
    SVLocusSetOptions sopt;
    sopt.minMergeEdgeObservations = 2;
    SVLocusSet set1(sopt);
    set1.merge(locus1);
    const SVLocusSet& cset1(set1);

    set1.cleanRegion(GenomeInterval(1, 5, 15));

    TestSVLocusSetProperties(cset1, 2, 0, 0, 0);
  }
}

BOOST_AUTO_TEST_CASE(test_SVLocusSmallDelRegionClean)
{
  // regions picked up from deletions have counts on both sides

  SVLocus locus1;
  locusAddPair(locus1, 1, 10, 20, 1, 30, 40, true);

  {
    SVLocusSetOptions sopt;
    sopt.minMergeEdgeObservations = 2;
    SVLocusSet set1(sopt);
    set1.merge(locus1);
    const SVLocusSet& cset1(set1);

    set1.cleanRegion(GenomeInterval(1, 0, 70));

    TestSVLocusSetProperties(cset1, 2, 0, 0, 0);
  }

  {
    SVLocusSetOptions sopt;
    sopt.minMergeEdgeObservations = 2;
    SVLocusSet set1(sopt);
    set1.merge(locus1);
    const SVLocusSet& cset1(set1);

    set1.clean();

    TestSVLocusSetProperties(cset1, 1, 0, 0, 0);
  }
}

BOOST_AUTO_TEST_CASE(test_SVLocusCleanSelfEdge)
{
  BOOST_TEST_MESSAGE("SDS MANTA-700");

  // Finalize with minMergeEdgeObservations clears out single self edge locus.

  SVLocus locus1;
  locusAddPair(locus1, 1, 10, 20, 1, 10, 20, true, 0);

  locus1.mergeSelfOverlap();

  {
    SVLocusSetOptions sopt;
    sopt.minMergeEdgeObservations = 3;
    SVLocusSet set1(sopt);
    set1.merge(locus1);
    const SVLocusSet& cset1(set1);

    TestSVLocusSetProperties(cset1, 1, 1, 1, 1);

    set1.finalize();
    cset1.checkState(true, true);
    TestSVLocusSetProperties(cset1, 1, 0, 0, 0);
  }
}

#if 0
BOOST_AUTO_TEST_CASE( test_SVLocusTransitiveOverlap )
{
    // abstracted from real-data error case on MANTA-28
    //
    // what happens when there's a complex transitive overlap chain?

    SVLocus locus1;
    locusAddPair(locus1,1,25,32,1,25,32,true,3);
    SVLocus locus2a;
    locusAddPair(locus2a,1,8,12,1,14,27,true,1);
    SVLocus locus2;
    locusAddPair(locus2,1,11,14,1,18,22,true,1);
    SVLocus locus3;
    locusAddPair(locus3,1,11,16,1,18,20,true,2);

    locus1.mergeSelfOverlap();
    locus2a.mergeSelfOverlap();
    locus2.mergeSelfOverlap();
    locus3.mergeSelfOverlap();
    {
        SVLocusSet set1(3);
        set1.merge(locus1);
        set1.merge(locus2a);
        set1.merge(locus2);
        set1.merge(locus3);
        const SVLocusSet& cset1(set1);

        BOOST_REQUIRE_EQUAL(cset1.nonEmptySize(),1u);
        set1.finalize();
        cset1.checkState(true,true);
        BOOST_REQUIRE_EQUAL(cset1.nonEmptySize(),1u);
    }
}
#endif

BOOST_AUTO_TEST_CASE(test_SVLocusTransitiveOverlap2)
{
  BOOST_TEST_MESSAGE("SDS MANTA-700");

  // abstracted from real-data error case on MANTA-28
  //
  // what happens when there's a complex transitive overlap chain?
#if 0
    // original coordinates extracted from actual failure case:
    SVLocus locus1;
    locusAddPair(locus1,1,615,837,1,853,900,true,6);
    SVLocus locus2a;
    locusAddPair(locus2a,1,464,614,1,712,862,true,1);
    SVLocus locus2b;
    locusAddPair(locus2b,1,645,798,1,421,574,false,1);
    SVLocus locus2c;
    locusAddPair(locus2c,1,370,851,1,370,851,true,3);
    SVLocus locus3;
    locusAddPair(locus3,1,693,843,1,538,688,false,1);
#endif

  SVLocus locus1;
  locusAddPair(locus1, 1, 30, 40, 1, 50, 60, true, 6);
  SVLocus locus2a;
  locusAddPair(locus2a, 1, 10, 20, 1, 30, 60, true, 1);
  SVLocus locus2b;
  locusAddPair(locus2b, 1, 30, 40, 1, 10, 20, false, 1);
  SVLocus locus2c;
  locusAddPair(locus2c, 1, 10, 40, 1, 10, 40, true, 3);
  SVLocus locus3;
  locusAddPair(locus3, 1, 30, 40, 1, 10, 20, false, 1);

  locus1.mergeSelfOverlap();
  locus2a.mergeSelfOverlap();
  locus2b.mergeSelfOverlap();
  locus2c.mergeSelfOverlap();
  locus3.mergeSelfOverlap();
  {
    SVLocusSetOptions sopt;
    sopt.minMergeEdgeObservations = 6;
    SVLocusSet set1(sopt);
    set1.merge(locus1);
    set1.merge(locus2c);
    set1.merge(locus2a);
    set1.merge(locus2b);
    set1.merge(locus3);
    const SVLocusSet& cset1(set1);

    TestSVLocusSetProperties(cset1, 5, 1, 1, 1);
    set1.finalize();
    cset1.checkState(true, true);
    TestSVLocusSetProperties(cset1, 5, 1, 1, 1);
    BOOST_REQUIRE_EQUAL(cset1.getLocus(0).getNode(0).getInterval(), GenomeInterval(1, 10, 60));
    BOOST_REQUIRE_EQUAL(cset1.getLocus(0).getNode(0).outCount(), 12);
  }
}

BOOST_AUTO_TEST_CASE(test_SVLocusTransitiveOverlap3)
{
  // abstracted from real-data error case on MANTA-28
  //
  // what happens when there's a complex transitive overlap chain?

  SVLocus locus1;
  locusAddPair(locus1, 1, 40, 60, 1, 70, 80, true, 2);
  SVLocus locus2a;
  locusAddPair(locus2a, 1, 10, 40, 1, 50, 60, true, 1);
  SVLocus locus3;
  locusAddPair(locus3, 1, 10, 20, 1, 30, 60, false, 1);

  locus1.mergeSelfOverlap();
  locus2a.mergeSelfOverlap();
  locus3.mergeSelfOverlap();
  {
    SVLocusSetOptions sopt;
    sopt.minMergeEdgeObservations = 2;
    SVLocusSet set1(sopt);
    set1.merge(locus1);
    set1.merge(locus2a);
    set1.merge(locus3);
    const SVLocusSet& cset1(set1);

    TestSVLocusSetProperties(cset1, 3, 1, 2, 3);
    set1.finalize();
    cset1.checkState(true, true);
    TestSVLocusSetProperties(cset1, 3, 1, 2, 3);
    BOOST_REQUIRE_EQUAL(cset1.getLocus(0).getNode(0).getInterval(), GenomeInterval(1, 10, 60));
    BOOST_REQUIRE_EQUAL(cset1.getLocus(0).getNode(1).getInterval(), GenomeInterval(1, 70, 80));
    BOOST_REQUIRE_EQUAL(cset1.getLocus(0).getNode(0).outCount(), 4);
    BOOST_REQUIRE_EQUAL(cset1.getLocus(0).getNode(1).outCount(), 2);
  }
}

// replicate at least one part of MANTA-257 in minimal form:
BOOST_AUTO_TEST_CASE(test_SVLocusSet_MANTA257_min1)
{
  SVLocus locus1;
  {
    const NodeIndexType nodePtr1 = locus1.addNode(GenomeInterval(0, 10, 20));
    const NodeIndexType nodePtr2 = locus1.addNode(GenomeInterval(1, 60, 80));
    const NodeIndexType nodePtr3 = locus1.addNode(GenomeInterval(1, 20, 50));

    locus1.linkNodes(nodePtr1, nodePtr2);
    locus1.linkNodes(nodePtr1, nodePtr3);
  }

  SVLocus locus2;
  {
    const NodeIndexType nodePtr1 = locus2.addNode(GenomeInterval(1, 10, 30));
    const NodeIndexType nodePtr2 = locus2.addNode(GenomeInterval(0, 10, 20));
    const NodeIndexType nodePtr3 = locus2.addNode(GenomeInterval(1, 40, 70));

    locus2.linkNodes(nodePtr1, nodePtr2);
    locus2.linkNodes(nodePtr3, nodePtr1);
  }

  SVLocusSetOptions sopt;
  sopt.minMergeEdgeObservations = 1;
  SVLocusSet set1(sopt);
  set1.merge(locus1);
  set1.merge(locus2);
  const SVLocusSet& cset1(set1);

  set1.finalize();
  cset1.checkState(true, true);
  TestSVLocusSetProperties(cset1, 2, 1, 2, 3);
  BOOST_REQUIRE_EQUAL(cset1.getLocus(0).getNode(0).getInterval(), GenomeInterval(0, 10, 20));
  BOOST_REQUIRE_EQUAL(cset1.getLocus(0).getNode(1).getInterval(), GenomeInterval(1, 10, 80));
  BOOST_REQUIRE_EQUAL(cset1.getLocus(0).getNode(0).outCount(), 2);
  BOOST_REQUIRE_EQUAL(cset1.getLocus(0).getNode(1).outCount(), 2);
}

// replicate the MANTA-257 bug in slightly reduced form:
BOOST_AUTO_TEST_CASE(test_SVLocusSet_MANTA257_simplified)
{
  SVLocus locus1;
  {
    const NodeIndexType nodePtr1 = locus1.addNode(GenomeInterval(0, 10, 40));
    const NodeIndexType nodePtr2 = locus1.addNode(GenomeInterval(1, 60, 100));
    const NodeIndexType nodePtr3 = locus1.addNode(GenomeInterval(1, 20, 50));

    locus1.linkNodes(nodePtr1, nodePtr2, 1, 0);
    locus1.linkNodes(nodePtr1, nodePtr3, 1, 0);
  }

  SVLocus locus2;
  {
    const NodeIndexType nodePtr1 = locus2.addNode(GenomeInterval(1, 10, 30));
    const NodeIndexType nodePtr2 = locus2.addNode(GenomeInterval(0, 20, 30));
    const NodeIndexType nodePtr3 = locus2.addNode(GenomeInterval(1, 80, 90));
    const NodeIndexType nodePtr4 = locus2.addNode(GenomeInterval(1, 40, 70));

    locus2.linkNodes(nodePtr1, nodePtr2, 1, 0);
    locus2.linkNodes(nodePtr1, nodePtr3, 1, 0);
    locus2.linkNodes(nodePtr4, nodePtr1, 1, 0);
  }

  SVLocusSetOptions sopt;
  sopt.minMergeEdgeObservations = 1;
  SVLocusSet set1(sopt);
  set1.merge(locus1);
  set1.merge(locus2);
  const SVLocusSet& cset1(set1);

  set1.finalize();
  cset1.checkState(true, true);
  TestSVLocusSetProperties(cset1, 2, 1, 2, 3);
  BOOST_REQUIRE_EQUAL(cset1.getLocus(0).getNode(0).getInterval(), GenomeInterval(0, 10, 40));
  BOOST_REQUIRE_EQUAL(cset1.getLocus(0).getNode(1).getInterval(), GenomeInterval(1, 10, 100));
  BOOST_REQUIRE_EQUAL(cset1.getLocus(0).getNode(0).outCount(), 2);
  BOOST_REQUIRE_EQUAL(cset1.getLocus(0).getNode(1).outCount(), 3);
}

// replicate the MANTA-257 bug in minimally reduced form:
BOOST_AUTO_TEST_CASE(test_SVLocusSet_MANTA257)
{
  SVLocus locus1;
  {
    const NodeIndexType nodePtr1 = locus1.addNode(GenomeInterval(0, 2255650, 2256356));
    const NodeIndexType nodePtr2 = locus1.addNode(GenomeInterval(1, 776, 1618));
    const NodeIndexType nodePtr3 = locus1.addNode(GenomeInterval(1, -298, 488));

    locus1.linkNodes(nodePtr1, nodePtr2, 51, 0);
    locus1.linkNodes(nodePtr1, nodePtr3, 78, 0);
  }

  SVLocus locus2;
  {
    const NodeIndexType nodePtr1 = locus2.addNode(GenomeInterval(1, -309, 265));
    const NodeIndexType nodePtr2 = locus2.addNode(GenomeInterval(0, 2255700, 2256245));
    const NodeIndexType nodePtr3 = locus2.addNode(GenomeInterval(1, 1018, 1595));
    const NodeIndexType nodePtr4 = locus2.addNode(GenomeInterval(1, 412, 904));

    locus2.linkNodes(nodePtr1, nodePtr2, 21, 0);
    locus2.linkNodes(nodePtr1, nodePtr3, 9, 3);
    locus2.linkNodes(nodePtr4, nodePtr1, 12, 0);
  }

  SVLocusSetOptions sopt;
  sopt.minMergeEdgeObservations = 9;
  SVLocusSet set1(sopt);
  set1.merge(locus1);
  set1.merge(locus2);
  const SVLocusSet& cset1(set1);

  set1.finalize();
  cset1.checkState(true, true);
  TestSVLocusSetProperties(cset1, 2, 1, 2, 3);
  BOOST_REQUIRE_EQUAL(cset1.getLocus(0).getNode(0).getInterval(), GenomeInterval(0, 2255650, 2256356));
  BOOST_REQUIRE_EQUAL(cset1.getLocus(0).getNode(1).getInterval(), GenomeInterval(1, -309, 1618));
  BOOST_REQUIRE_EQUAL(cset1.getLocus(0).getNode(0).outCount(), 129);
  BOOST_REQUIRE_EQUAL(cset1.getLocus(0).getNode(1).outCount(), 42);
}

BOOST_AUTO_TEST_CASE(test_SVLocusSet_Save_Load)
{
  BOOST_TEST_MESSAGE("SDS MANTA-701");
  BOOST_TEST_MESSAGE("SDS MANTA-702");

  // construct a simple two-node locus
  SVLocus locus1;
  locusAddPair(locus1, 1, 10, 20, 2, 30, 40);

  SVLocus locus2;
  locusAddPair(locus2, 3, 60, 80, 4, 90, 100);

  SVLocus locus3;
  locusAddPair(locus3, 1, 10, 20, 2, 30, 40);

  SVLocusSetOptions sopt;
  sopt.minMergeEdgeObservations = 1;

  SVLocusSet set1(sopt);

  set1.merge(locus1);
  set1.merge(locus3);
  set1.merge(locus2);
  set1.checkState(true, true);

  const auto set1_copy_ptr(getSerializedSVLocusSetCopy(set1));

  SVLocusSet& cset1(set1);
  SVLocusSet& cset1_copy(*set1_copy_ptr);

  // Test if SVLocusSet2 == SVLocusSet1
  cset1_copy.checkState(true, true);
  TestSVLocusSetProperties(cset1, 2, 2, 4, 4);
  TestSVLocusSetProperties(cset1_copy, 2, 2, 4, 4);
}

BOOST_AUTO_TEST_CASE(test_SVLocusSet_DumpLoci)
{
  // construct a simple two-node locus
  SVLocus locus1;
  locusAddPair(locus1, 1, 10, 20, 2, 30, 40);

  SVLocus locus2;
  locusAddPair(locus2, 1, 10, 20, 2, 30, 40);

  SVLocusSetOptions sopt;
  sopt.minMergeEdgeObservations = 2;
  SVLocusSet set1(sopt);

  set1.merge(locus1);
  set1.merge(locus2);
  set1.checkState(true, true);

  // print stats regarding the SVLocusSet.
  std::string  test1 = "SVLocusSetTest1.txt";
  std::filebuf fileBuffer;
  fileBuffer.open(test1, std::ios::out);
  std::ostream outStream1(&fileBuffer);
  set1.dumpStats(outStream1);
  fileBuffer.close();
  BOOST_REQUIRE(FindStringInFile(test1, "GraphBuildTime"));
  BOOST_REQUIRE(FindStringInFile(test1, "directedEdges"));
  std::remove(test1.c_str());

  // Print stats regarding the loci in the SVLocusSet.
  std::string test2 = "SVLocusSetTest2.txt";
  fileBuffer.open(test2, std::ios::out);
  std::ostream outStream2(&fileBuffer);
  set1.dumpLocusStats(outStream2);
  fileBuffer.close();
  BOOST_REQUIRE(FindStringInFile(test2, "locusIndex"));
  BOOST_REQUIRE(FindStringInFile(test2, "maxNodeObsCount"));
  std::remove(test2.c_str());

  // Print all of the loci in the SVLocusset.
  std::string test3 = "SVLocusSetTest3.txt";
  fileBuffer.open(test3, std::ios::out);
  std::ostream outStream3(&fileBuffer);
  set1.dump(outStream3);
  fileBuffer.close();
  BOOST_REQUIRE(FindStringInFile(test3, "LOCUSSET_START"));
  BOOST_REQUIRE(FindStringInFile(
      test3,
      "NodeIndex: 0 LocusNode: GenomeInterval: 1:[10,20) n_edges: 1 out_count: 2 in_count: 0 evidence: [10,20)"));
  BOOST_REQUIRE(FindStringInFile(
      test3,
      "NodeIndex: 1 LocusNode: GenomeInterval: 2:[30,40) n_edges: 1 out_count: 0 in_count: 2 evidence: [30,40)"));
  BOOST_REQUIRE(FindStringInFile(test3, "LOCUSSET_END"));
  std::remove(test3.c_str());

  // Print a region of loci. Does not have section header/footer
  std::string test4 = "SVLocusSetTest4.txt";
  fileBuffer.open(test4, std::ios::out);
  std::ostream         outStream4(&fileBuffer);
  const GenomeInterval dumpRegion(1, 10, 30);
  set1.dumpRegion(outStream4, dumpRegion);
  fileBuffer.close();
  BOOST_REQUIRE(!FindStringInFile(test4, "LOCUSSET_START"));
  BOOST_REQUIRE(!FindStringInFile(
      test4,
      "NodeIndex: 0 LocusNode: GenomeInterval: 1:[10,20) n_edges: 1 out_count: 2 in_count: 0 evidence: [10,20)"));
  BOOST_REQUIRE(FindStringInFile(
      test4, "LocusNode: GenomeInterval: 1:[10,20) n_edges: 1 out_count: 2 evidence: [10,20)"));
  BOOST_REQUIRE(!FindStringInFile(test4, "LOCUSSET_END"));
  std::remove(test4.c_str());
}

BOOST_AUTO_TEST_SUITE_END()
