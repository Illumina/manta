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

#include "svgraph/GenomeInterval.hpp"

// Test public functionality of the GenomeInterval
BOOST_AUTO_TEST_SUITE(test_GenomeInterval)

BOOST_AUTO_TEST_CASE(test_GenomeInterval_Clear)
{
  GenomeInterval TestInterval1(1, 100, 200);
  BOOST_REQUIRE_EQUAL(TestInterval1.tid, 1);
  BOOST_REQUIRE_EQUAL(TestInterval1.range.begin_pos(), 100);
  BOOST_REQUIRE_EQUAL(TestInterval1.range.end_pos(), 200);
  TestInterval1.clear();
  BOOST_REQUIRE_EQUAL(TestInterval1.tid, 0);
  BOOST_REQUIRE_EQUAL(TestInterval1.range.begin_pos(), 0);
  BOOST_REQUIRE_EQUAL(TestInterval1.range.end_pos(), 0);
}

BOOST_AUTO_TEST_CASE(test_GenomeInterval_LessThanOperator)
{
  // < operator tests
  BOOST_REQUIRE_EQUAL(GenomeInterval(1, 10, 40) < GenomeInterval(1, 10, 40), false);

  // if id less than, then true, regardless of range
  BOOST_REQUIRE_EQUAL(GenomeInterval(1, 10, 40) < GenomeInterval(2, 10, 40), true);
  BOOST_REQUIRE_EQUAL(GenomeInterval(1, 10, 40) < GenomeInterval(2, 41, 50), true);
  BOOST_REQUIRE_EQUAL(GenomeInterval(1, 10, 40) < GenomeInterval(2, 5, 10), true);

  // if id is the same, test for range
  BOOST_REQUIRE_EQUAL(GenomeInterval(1, 10, 40) < GenomeInterval(1, 5, 9), false);
  BOOST_REQUIRE_EQUAL(GenomeInterval(1, 10, 40) < GenomeInterval(1, 41, 50), true);
  BOOST_REQUIRE_EQUAL(GenomeInterval(1, 10, 40) < GenomeInterval(1, 5, 10), false);
  BOOST_REQUIRE_EQUAL(GenomeInterval(1, 10, 40) < GenomeInterval(1, 20, 40), true);
  BOOST_REQUIRE_EQUAL(GenomeInterval(1, 10, 40) < GenomeInterval(1, 9, 40), false);
  BOOST_REQUIRE_EQUAL(GenomeInterval(1, 10, 40) < GenomeInterval(1, 11, 40), true);
  BOOST_REQUIRE_EQUAL(GenomeInterval(1, 10, 40) < GenomeInterval(1, 10, 39), false);
  BOOST_REQUIRE_EQUAL(GenomeInterval(1, 10, 40) < GenomeInterval(1, 10, 41), true);
  BOOST_REQUIRE_EQUAL(GenomeInterval(1, 10, 40) < GenomeInterval(1, 11, 39), true);
  BOOST_REQUIRE_EQUAL(GenomeInterval(1, 10, 40) < GenomeInterval(1, 9, 41), false);
}

BOOST_AUTO_TEST_CASE(test_GenomeInterval_EqualOperator)
{
  // == operator tests, all values must be identical
  BOOST_REQUIRE_EQUAL(GenomeInterval(1, 10, 40) == GenomeInterval(1, 10, 40), true);
  BOOST_REQUIRE_EQUAL(GenomeInterval(1, 10, 40) == GenomeInterval(2, 10, 40), false);
  BOOST_REQUIRE_EQUAL(GenomeInterval(1, 10, 40) == GenomeInterval(1, 20, 40), false);
  BOOST_REQUIRE_EQUAL(GenomeInterval(1, 10, 40) == GenomeInterval(1, 10, 39), false);
}

BOOST_AUTO_TEST_CASE(test_GenomeInterval_Intersetions)
{
  GenomeInterval genomeBaseInterval(1, 10, 40);
  GenomeInterval genomeIntervalDifferentId(2, 10, 40);
  GenomeInterval genomeIntervalOverlapping(1, 1, 100);
  GenomeInterval genomeIntervalBefore(1, 1, 5);
  GenomeInterval genomeIntervalAfter(1, 90, 100);
  GenomeInterval genomeInterval1BPOverlapBeginning(1, 1, 10);
  GenomeInterval genomeInterval1BPOverlapEnding(1, 40, 100);
  GenomeInterval genomeInterval2BPOverlapBeginning(1, 1, 11);
  GenomeInterval genomeInterval2BPOverlapEnding(1, 39, 100);
  GenomeInterval genomeIntervalInside(1, 11, 39);

  // Test Genome Interval Intersection
  BOOST_REQUIRE_EQUAL(genomeBaseInterval.isIntersect(genomeBaseInterval), true);
  BOOST_REQUIRE_EQUAL(genomeBaseInterval.isIntersect(genomeIntervalDifferentId), false);
  BOOST_REQUIRE_EQUAL(genomeBaseInterval.isIntersect(genomeIntervalOverlapping), true);
  BOOST_REQUIRE_EQUAL(genomeBaseInterval.isIntersect(genomeIntervalBefore), false);
  BOOST_REQUIRE_EQUAL(genomeBaseInterval.isIntersect(genomeIntervalAfter), false);

  BOOST_REQUIRE_EQUAL(genomeBaseInterval.isIntersect(genomeInterval1BPOverlapBeginning), false);
  BOOST_REQUIRE_EQUAL(genomeBaseInterval.isIntersect(genomeInterval2BPOverlapBeginning), true);

  BOOST_REQUIRE_EQUAL(genomeBaseInterval.isIntersect(genomeInterval1BPOverlapEnding), false);
  BOOST_REQUIRE_EQUAL(genomeBaseInterval.isIntersect(genomeInterval2BPOverlapEnding), true);

  BOOST_REQUIRE_EQUAL(genomeBaseInterval.isIntersect(genomeIntervalInside), true);
}

BOOST_AUTO_TEST_CASE(test_GenomeInterval_Sorting)
{
  // test that GenomeInterval sorting:
  std::vector<GenomeInterval> test;

  test.push_back(GenomeInterval(1, 15, 19));
  test.push_back(GenomeInterval(1, 15, 22));
  test.push_back(GenomeInterval(1, 10, 20));
  test.push_back(GenomeInterval(2, 5, 10));
  test.push_back(GenomeInterval(2, 8, 10));

  std::sort(test.begin(), test.end());

  BOOST_REQUIRE_EQUAL(test[0], GenomeInterval(1, 10, 20));
  BOOST_REQUIRE_EQUAL(test[2], GenomeInterval(1, 15, 22));
  BOOST_REQUIRE_EQUAL(test[4], GenomeInterval(2, 8, 10));
}

#if 0
BOOST_AUTO_TEST_CASE( test_GenomeInterval_OStreamOperator )
{
    // TODO Need a mock framework to do this cleanly by mocking the streambuffer
}
#endif

BOOST_AUTO_TEST_SUITE_END()
