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

#include <fstream>
#include <thread>

#include "EdgeRuntimeTracker.hpp"
#include "boost/filesystem.hpp"
#include "boost/test/unit_test.hpp"
#include "test/testFileMakers.hpp"

/// TestEdgeRuntimeTracker is a friend of EdgeRuntimeTracker. So that can access private
/// members of EdgeRuntimeTracker.
struct TestEdgeRuntimeTracker {
  unsigned getCandidate(EdgeRuntimeTracker& tracker) { return tracker._candidateCount; }

  unsigned getComplexCandidate(EdgeRuntimeTracker& tracker) { return tracker._complexCandidateCount; }

  unsigned getAssembledCandidate(EdgeRuntimeTracker& tracker) { return tracker._assembledCandidateCount; }

  unsigned getAssembledComplexCandidate(EdgeRuntimeTracker& tracker)
  {
    return tracker._assembledComplexCandidateCount;
  }
};

BOOST_AUTO_TEST_SUITE(EdgeRuntimeTracket_test_suite)

// Following statistics are verifed from EdgeTracker file
// 1. Edge info (locusIndex:nodeindex1:nodeIndex2)
// 2. candidate count
// 3. complex candidate count
// 4. Assembled candidate count
// 5. Complex assembled candidate count
BOOST_AUTO_TEST_CASE(test_tracker)
{
  TestFilenameMaker  filenameMaker;
  EdgeRuntimeTracker tracker(filenameMaker.getFilename());
  EdgeInfo           info;
  info.nodeIndex1 = 1;
  info.nodeIndex2 = 2;
  tracker.start();
  // Adding 2 assembled candidates, 1 complex assembled candidate,
  // 3 candidates, 1 complex candidates.
  tracker.addAssembledCandidate(false);
  tracker.addAssembledCandidate(false);
  tracker.addAssembledCandidate(true);
  tracker.addCandidate(false);
  tracker.addCandidate(false);
  tracker.addCandidate(false);
  tracker.addCandidate(true);
  // Adding 1 sec delay as for writing at least 0.5 sec interval is required.
  std::this_thread::sleep_for(std::chrono::seconds(1));
  tracker.stop(info);
  TestEdgeRuntimeTracker testEdgeRuntimeTracker;

  // Verify the result
  BOOST_REQUIRE_EQUAL(testEdgeRuntimeTracker.getCandidate(tracker), 3);
  BOOST_REQUIRE_EQUAL(testEdgeRuntimeTracker.getComplexCandidate(tracker), 1);
  BOOST_REQUIRE_EQUAL(testEdgeRuntimeTracker.getAssembledCandidate(tracker), 2);
  BOOST_REQUIRE_EQUAL(testEdgeRuntimeTracker.getAssembledComplexCandidate(tracker), 1);
}

BOOST_AUTO_TEST_SUITE_END()
