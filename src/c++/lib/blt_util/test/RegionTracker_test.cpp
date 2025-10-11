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

#include "RegionTracker.hpp"

BOOST_AUTO_TEST_SUITE(test_RegionTracker)

BOOST_AUTO_TEST_CASE(test_RegionTrackerSimple)
{
  // Simplest test
  RegionTracker rt;

  rt.addRegion(known_pos_range2(0, 1));
  BOOST_REQUIRE(rt.isIntersectRegion(0));
  BOOST_REQUIRE(!rt.isIntersectRegion(1));
}

BOOST_AUTO_TEST_CASE(test_RegionTrackerPosIntersect)
{
  // region overlap tests
  {
    RegionTracker rt;

    rt.addRegion(known_pos_range2(5, 10));
    rt.addRegion(known_pos_range2(2, 3));
    BOOST_REQUIRE_EQUAL(rt.size(), 2);
    BOOST_REQUIRE(rt.isIntersectRegion(2));
    BOOST_REQUIRE(!rt.isIntersectRegion(3));
    BOOST_REQUIRE(!rt.isIntersectRegion(4));
    BOOST_REQUIRE(rt.isIntersectRegion(5));

    rt.addRegion(known_pos_range2(3, 7));
    BOOST_REQUIRE_EQUAL(rt.size(), 1);
    BOOST_REQUIRE(rt.isIntersectRegion(4));
  }
  {
    RegionTracker rt;
    rt.addRegion(known_pos_range2(5, 10));
    rt.addRegion(known_pos_range2(2, 3));
    rt.addRegion(known_pos_range2(2, 5));
    BOOST_REQUIRE_EQUAL(rt.size(), 1);
    BOOST_REQUIRE(rt.isIntersectRegion(4));
  }

  {
    RegionTracker rt;
    rt.addRegion(known_pos_range2(5, 10));
    rt.addRegion(known_pos_range2(2, 3));
    rt.addRegion(known_pos_range2(2, 4));
    BOOST_REQUIRE_EQUAL(rt.size(), 2);
    BOOST_REQUIRE(rt.isIntersectRegion(3));
  }

  {
    RegionTracker rt;
    rt.addRegion(known_pos_range2(5, 10));
    rt.addRegion(known_pos_range2(2, 3));
    rt.addRegion(known_pos_range2(4, 5));
    BOOST_REQUIRE_EQUAL(rt.size(), 2);
    BOOST_REQUIRE(rt.isIntersectRegion(4));
  }

  {
    RegionTracker rt;
    rt.addRegion(known_pos_range2(1, 10));
    rt.addRegion(known_pos_range2(4, 5));
    BOOST_REQUIRE_EQUAL(rt.size(), 1);
    BOOST_REQUIRE(rt.isIntersectRegion(4));
  }

  {
    RegionTracker rt;
    rt.addRegion(known_pos_range2(4, 5));
    rt.addRegion(known_pos_range2(1, 10));
    BOOST_REQUIRE_EQUAL(rt.size(), 1);
    BOOST_REQUIRE(rt.isIntersectRegion(4));
  }
}

BOOST_AUTO_TEST_CASE(test_RegionTrackerRegionIntersect)
{
  // region overlap tests
  {
    RegionTracker rt;

    rt.addRegion(known_pos_range2(5, 10));
    rt.addRegion(known_pos_range2(2, 3));
    BOOST_REQUIRE_EQUAL(rt.size(), 2);
    BOOST_REQUIRE(rt.isIntersectRegion(known_pos_range2(2, 10)));
    BOOST_REQUIRE(!rt.isIntersectRegion(known_pos_range2(3, 4)));
    BOOST_REQUIRE(!rt.isIntersectRegion(known_pos_range2(4, 4)));
    BOOST_REQUIRE(rt.isIntersectRegion(known_pos_range2(5, 11)));
  }
}

BOOST_AUTO_TEST_CASE(test_RegionTrackerRemoveToPos)
{
  // region remove tests
  RegionTracker rt;

  rt.addRegion(known_pos_range2(5, 10));
  rt.addRegion(known_pos_range2(2, 3));
  rt.addRegion(known_pos_range2(14, 15));
  rt.addRegion(known_pos_range2(24, 25));
  BOOST_REQUIRE_EQUAL(rt.size(), 4);

  rt.removeToPos(2);
  BOOST_REQUIRE_EQUAL(rt.size(), 3);
  rt.removeToPos(2);
  BOOST_REQUIRE_EQUAL(rt.size(), 3);
  rt.removeToPos(6);
  BOOST_REQUIRE_EQUAL(rt.size(), 3);
  rt.removeToPos(16);
  BOOST_REQUIRE_EQUAL(rt.size(), 1);
}

BOOST_AUTO_TEST_CASE(test_RegionTrackerSubset)
{
  // region remove tests
  RegionTracker rt;

  rt.addRegion(known_pos_range2(5, 10));
  rt.addRegion(known_pos_range2(2, 3));
  rt.addRegion(known_pos_range2(14, 15));
  rt.addRegion(known_pos_range2(24, 25));
  BOOST_REQUIRE_EQUAL(rt.size(), 4);

  BOOST_REQUIRE(rt.isSubsetOfRegion(known_pos_range2(5, 10)));
  BOOST_REQUIRE(rt.isSubsetOfRegion(known_pos_range2(6, 7)));
  BOOST_REQUIRE(!rt.isSubsetOfRegion(known_pos_range2(4, 10)));
  BOOST_REQUIRE(!rt.isSubsetOfRegion(known_pos_range2(5, 11)));
  BOOST_REQUIRE(!rt.isSubsetOfRegion(known_pos_range2(0, 1)));
  BOOST_REQUIRE(!rt.isSubsetOfRegion(known_pos_range2(30, 31)));
}

BOOST_AUTO_TEST_CASE(test_RegionTrackerMerge)
{
  RegionTracker rt1, rt2;

  rt1.addRegion(known_pos_range2(5, 10));
  rt1.addRegion(known_pos_range2(14, 15));

  rt2.addRegion(known_pos_range2(2, 3));
  rt2.addRegion(known_pos_range2(24, 25));

  rt1.merge(rt2);

  BOOST_REQUIRE_EQUAL(rt1.size(), 4);
  BOOST_REQUIRE_EQUAL(rt2.size(), 2);

  BOOST_REQUIRE(rt1.isIntersectRegion(24));
  BOOST_REQUIRE(rt1.isSubsetOfRegion(known_pos_range2(5, 10)));
  BOOST_REQUIRE(rt1.isSubsetOfRegion(known_pos_range2(6, 7)));
  BOOST_REQUIRE(!rt1.isSubsetOfRegion(known_pos_range2(4, 10)));
  BOOST_REQUIRE(!rt1.isSubsetOfRegion(known_pos_range2(5, 11)));
  BOOST_REQUIRE(!rt1.isSubsetOfRegion(known_pos_range2(0, 1)));
  BOOST_REQUIRE(!rt1.isSubsetOfRegion(known_pos_range2(30, 31)));
}

BOOST_AUTO_TEST_CASE(test_RegionPayloadTracker)
{
  // Simplest test
  RegionPayloadTracker<int> rt;

  rt.addRegion(known_pos_range2(0, 1), 25);
  auto val = rt.isIntersectRegion(0);
  BOOST_REQUIRE(val);
  BOOST_REQUIRE_EQUAL(*val, 25);
  BOOST_REQUIRE(!rt.isIntersectRegion(1));
}

BOOST_AUTO_TEST_CASE(test_RegionPayloadTracker2)
{
  // region overlap tests
  {
    RegionPayloadTracker<int> rt;

    BOOST_REQUIRE(rt.addRegion(known_pos_range2(5, 10), 5));
    BOOST_REQUIRE(rt.addRegion(known_pos_range2(2, 3), 5));
    BOOST_REQUIRE_EQUAL(rt.size(), 2);
    BOOST_REQUIRE(rt.isIntersectRegion(2));
    BOOST_REQUIRE(!rt.isIntersectRegion(3));
    BOOST_REQUIRE(!rt.isIntersectRegion(4));
    BOOST_REQUIRE(rt.isIntersectRegion(5));

    BOOST_REQUIRE(rt.addRegion(known_pos_range2(3, 7), 5));
    BOOST_REQUIRE_EQUAL(rt.size(), 1);
    BOOST_REQUIRE(rt.isIntersectRegion(4));
  }
  {
    RegionPayloadTracker<int> rt;
    BOOST_REQUIRE(rt.addRegion(known_pos_range2(5, 10), 5));
    BOOST_REQUIRE(rt.addRegion(known_pos_range2(2, 3), 5));
    BOOST_REQUIRE(rt.addRegion(known_pos_range2(2, 5), 5));
    BOOST_REQUIRE_EQUAL(rt.size(), 1);
    BOOST_REQUIRE(rt.isIntersectRegion(4));
  }

  {
    RegionPayloadTracker<int> rt;
    rt.addRegion(known_pos_range2(5, 10), 5);
    rt.addRegion(known_pos_range2(2, 3), 5);
    rt.addRegion(known_pos_range2(2, 4), 5);
    BOOST_REQUIRE_EQUAL(rt.size(), 2);
    BOOST_REQUIRE(rt.isIntersectRegion(3));
  }

  {
    RegionPayloadTracker<int> rt;
    rt.addRegion(known_pos_range2(5, 10), 5);
    rt.addRegion(known_pos_range2(2, 3), 5);
    rt.addRegion(known_pos_range2(4, 5), 5);
    BOOST_REQUIRE_EQUAL(rt.size(), 2);
    BOOST_REQUIRE(rt.isIntersectRegion(4));
  }

  {
    RegionPayloadTracker<int> rt;
    rt.addRegion(known_pos_range2(1, 10), 5);
    rt.addRegion(known_pos_range2(4, 5), 5);
    BOOST_REQUIRE_EQUAL(rt.size(), 1);
    BOOST_REQUIRE(rt.isIntersectRegion(4));
  }

  {
    RegionPayloadTracker<int> rt;
    rt.addRegion(known_pos_range2(4, 5), 5);
    rt.addRegion(known_pos_range2(1, 10), 5);
    BOOST_REQUIRE_EQUAL(rt.size(), 1);
    BOOST_REQUIRE(rt.isIntersectRegion(4));
  }
}

BOOST_AUTO_TEST_CASE(test_RegionPayloadTracker3)
{
  // region remove tests
  RegionPayloadTracker<int> rt;

  rt.addRegion(known_pos_range2(5, 10), 5);
  rt.addRegion(known_pos_range2(2, 3), 5);
  rt.addRegion(known_pos_range2(14, 15), 5);
  rt.addRegion(known_pos_range2(24, 25), 5);
  BOOST_REQUIRE_EQUAL(rt.size(), 4);

  rt.removeToPos(2);
  BOOST_REQUIRE_EQUAL(rt.size(), 3);
  rt.removeToPos(2);
  BOOST_REQUIRE_EQUAL(rt.size(), 3);
  rt.removeToPos(6);
  BOOST_REQUIRE_EQUAL(rt.size(), 3);
  rt.removeToPos(16);
  BOOST_REQUIRE_EQUAL(rt.size(), 1);
}

BOOST_AUTO_TEST_CASE(test_RegionPayloadTracker4)
{
  // payload conflict tests
  {
    RegionPayloadTracker<int> rt;

    BOOST_REQUIRE(rt.addRegion(known_pos_range2(5, 10), 5));
    BOOST_REQUIRE(!rt.addRegion(known_pos_range2(8, 14), 4));
    BOOST_REQUIRE_EQUAL(rt.size(), 1);
    BOOST_REQUIRE(rt.addRegion(known_pos_range2(10, 14), 4));
    BOOST_REQUIRE(rt.addRegion(known_pos_range2(3, 5), 3));

    BOOST_REQUIRE(rt.isIntersectRegion(3));
    BOOST_REQUIRE_EQUAL(*rt.isIntersectRegion(3), 3);
  }
}

BOOST_AUTO_TEST_CASE(test_PayloadRegionTrackerSubset)
{
  // region remove tests
  RegionPayloadTracker<int> rt;

  rt.addRegion(known_pos_range2(5, 10), 5);
  rt.addRegion(known_pos_range2(2, 3), 4);
  rt.addRegion(known_pos_range2(14, 15), -1);
  rt.addRegion(known_pos_range2(24, 25), 0);
  BOOST_REQUIRE_EQUAL(rt.size(), 4);

  const auto test1(rt.isSubsetOfRegion(known_pos_range2(5, 10)));
  BOOST_REQUIRE(test1);
  BOOST_REQUIRE_EQUAL(*test1, 5);

  const auto test2(rt.isSubsetOfRegion(known_pos_range2(6, 7)));
  BOOST_REQUIRE(test2);
  BOOST_REQUIRE_EQUAL(*test2, 5);
  BOOST_REQUIRE(!rt.isSubsetOfRegion(known_pos_range2(4, 10)));
  BOOST_REQUIRE(!rt.isSubsetOfRegion(known_pos_range2(5, 11)));
  BOOST_REQUIRE(!rt.isSubsetOfRegion(known_pos_range2(0, 1)));
  BOOST_REQUIRE(!rt.isSubsetOfRegion(known_pos_range2(30, 31)));
}

BOOST_AUTO_TEST_SUITE_END()
