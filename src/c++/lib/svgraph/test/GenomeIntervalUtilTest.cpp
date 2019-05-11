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
