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

#include "svgraph/GenomeIntervalUtil.hpp"
#include "test/testAlignmentDataUtil.hpp"

BOOST_AUTO_TEST_SUITE(GenomeIntervalUtil_test_suite)

BOOST_AUTO_TEST_CASE(test_IntervalCompressor)
{
  // test that GenomeInterval sorting follows expect:
  std::vector<GenomeInterval> test;

  test.push_back(GenomeInterval(1, 15, 19));
  test.push_back(GenomeInterval(2, 5, 10));
  test.push_back(GenomeInterval(1, 10, 20));
  test.push_back(GenomeInterval(1, 24, 50));
  test.push_back(GenomeInterval(2, 8, 10));
  test.push_back(GenomeInterval(1, 15, 22));

  const std::vector<unsigned> indexMap = intervalCompressor(test);

  BOOST_REQUIRE_EQUAL(test.size(), 3u);
  BOOST_REQUIRE_EQUAL(test[0], GenomeInterval(1, 10, 22));
  BOOST_REQUIRE_EQUAL(test[1], GenomeInterval(2, 5, 10));

  BOOST_REQUIRE_EQUAL(indexMap.size(), 6u);
  BOOST_REQUIRE_EQUAL(indexMap[0], 0u);
  BOOST_REQUIRE_EQUAL(indexMap[5], 0u);
  BOOST_REQUIRE_EQUAL(indexMap[4], 1u);
}

BOOST_AUTO_TEST_CASE(test_convertSamtoolsRegionToGenomeInterval)
{
  BOOST_REQUIRE_EQUAL(
      GenomeInterval(0, 0, 500), convertSamtoolsRegionToGenomeInterval(buildTestBamHeader(), "chrFoo"));

  BOOST_REQUIRE_EQUAL(
      GenomeInterval(0, 99, 200),
      convertSamtoolsRegionToGenomeInterval(buildTestBamHeader(), "chrFoo:100-200"));
}

BOOST_AUTO_TEST_SUITE_END()
