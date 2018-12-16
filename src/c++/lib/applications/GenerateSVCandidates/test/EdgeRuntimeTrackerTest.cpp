//
// Manta - Structural Variant and Indel Caller
// Copyright (c) 2013-2018 Illumina, Inc.
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

#include "boost/test/unit_test.hpp"
#include "test/testFileMakers.hh"
#include "boost/filesystem.hpp"
#include "EdgeRuntimeTracker.hh"

/// TestEdgeRuntimeTracker is a friend of EdgeRuntimeTracker. So that can access private
/// members of EdgeRuntimeTracker.
struct TestEdgeRuntimeTracker
{
    void flushStream(EdgeRuntimeTracker& tracker)
    {
        tracker._osPtr->flush();
    }
};

BOOST_AUTO_TEST_SUITE( EdgeRuntimeTracket_test_suite )

// Following statistics are verifed from EdgeTracker file
// 1. Edge info (locusIndex:nodeindex1:nodeIndex2)
// 2. candidate count
// 3. complex candidate count
// 4. Assembled candidate count
// 5. Complex assembled candidate count
BOOST_AUTO_TEST_CASE( test_tracker)
{
    TestFilenameMaker filenameMaker;
    EdgeRuntimeTracker tracker(filenameMaker.getFilename());
    EdgeInfo info;
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
    testEdgeRuntimeTracker.flushStream(tracker);
    std::ifstream trackerFile(filenameMaker.getFilename());
    std::string edge;
    std::string edgeTime;
    unsigned candidate;
    unsigned complexCandidate;
    unsigned assembledCandidate;
    unsigned assembledComplexCandidate;

    // Verify the result
    trackerFile >> edge >> edgeTime >> candidate >> complexCandidate
                >> assembledCandidate >> assembledComplexCandidate;
    BOOST_REQUIRE_EQUAL(edge, "0:1:2");
    BOOST_REQUIRE_EQUAL(candidate, 3);
    BOOST_REQUIRE_EQUAL(complexCandidate, 1);
    BOOST_REQUIRE_EQUAL(assembledCandidate, 2);
    BOOST_REQUIRE_EQUAL(assembledComplexCandidate, 1);
}

BOOST_AUTO_TEST_SUITE_END()