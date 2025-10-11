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

#include "RegionSum.hpp"

BOOST_AUTO_TEST_SUITE(test_RegionSum)

BOOST_AUTO_TEST_CASE(RegionSum_test)
{
  RegionSum<unsigned> rs;
  rs.add(known_pos_range2(3, 7), 1u);
  rs.add(known_pos_range2(4, 5), 2u);

  BOOST_REQUIRE_EQUAL(rs.maxVal(), 3u);
}

BOOST_AUTO_TEST_SUITE_END()
