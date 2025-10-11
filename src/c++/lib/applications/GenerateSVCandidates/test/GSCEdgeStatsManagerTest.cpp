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

#include "GSCEdgeStatsManager.hpp"
#include "boost/filesystem.hpp"
#include "boost/test/unit_test.hpp"
#include "fstream"
#include "test/testFileMakers.hpp"

BOOST_AUTO_TEST_SUITE(GSCEdgeStatsManager_test_suite)

// Following statistics are verifed from EdgeStats file
// 1. Total number of junction assembly overlaps have been skipped
// 2. Total input edge count
// 3. Number of candidates
// 4. Number of complex candidates
// 5. Number of filtered spanning candidates
// 6. Number of junctions
// 7. Number of complex junctions
// 8. Number of assembly candidates
// 9. Number of spanning assembly candidates
BOOST_AUTO_TEST_CASE(test_GSCEdgeStatsManager)
{
  TestFilenameMaker   filenameMaker2;
  GSCEdgeStatsManager edgeStatsManager;
  EdgeRuntimeTracker  tracker(filenameMaker2.getFilename());
  EdgeInfo            edgeInfo;
  edgeInfo.nodeIndex1 = 1;
  edgeInfo.nodeIndex2 = 2;
  SVFinderStats finderStats;
  // Increment input edge count by 1 and increment candidate count by 5.
  edgeStatsManager.updateEdgeCandidates(edgeInfo, 5, finderStats);
  // Increment assembly candidates and spanning assembly candidates by 3.
  edgeStatsManager.updateAssemblyCount(edgeInfo, 3, true);
  // Increment assembly candidates and spanning assembly candidates by 3.
  // Increment totalJunctionAssemblyOverlapSkips by 1.
  edgeStatsManager.updateAssemblyCount(edgeInfo, 3, true, true);
  // Increment junction count by 10
  edgeStatsManager.updateJunctionCandidateCounts(edgeInfo, 10, false);
  // Increment total complex candidate by 2 and increment total Spanning Candidate Filter count by 4
  edgeStatsManager.updateMJFilter(edgeInfo, 2, 4);
  // Update the times
  edgeStatsManager.updateScoredEdgeTime(edgeInfo, tracker);
  tracker.stop(edgeInfo);

  // put edgeStats through serialize/deserialize cycle to test these functions as well:
  TestFilenameMaker filenameMaker1;
  {
    GSCEdgeStats edgeStats(edgeStatsManager.returnStats());
    edgeStats.save(filenameMaker1.getFilename().c_str());
  }

  GSCEdgeStats edgeStats;
  // loading edge stats from the file
  edgeStats.load(filenameMaker1.getFilename().c_str());
  // Check all the counts according to the description
  BOOST_REQUIRE_EQUAL(edgeStats.edgeData.remoteEdges.totalJunctionAssemblyOverlapSkips, 1);
  BOOST_REQUIRE_EQUAL(edgeStats.edgeData.remoteEdges.totalInputEdgeCount, 1);
  BOOST_REQUIRE_EQUAL(edgeStats.edgeData.remoteEdges.totalCandidateCount, 5);
  BOOST_REQUIRE_EQUAL(edgeStats.edgeData.remoteEdges.totalComplexCandidate, 2);
  BOOST_REQUIRE_EQUAL(edgeStats.edgeData.remoteEdges.totalSpanningCandidateFilter, 4);
  BOOST_REQUIRE_EQUAL(edgeStats.edgeData.remoteEdges.totalJunctionCount, 10);
  BOOST_REQUIRE_EQUAL(edgeStats.edgeData.remoteEdges.totalComplexJunctionCount, 0);
  BOOST_REQUIRE_EQUAL(edgeStats.edgeData.remoteEdges.totalAssemblyCandidates, 6);
  BOOST_REQUIRE_EQUAL(edgeStats.edgeData.remoteEdges.totalSpanningAssemblyCandidates, 6);
}

BOOST_AUTO_TEST_SUITE_END()
