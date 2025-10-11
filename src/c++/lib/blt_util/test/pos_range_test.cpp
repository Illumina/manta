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

#include "blt_util/pos_range.hpp"

BOOST_AUTO_TEST_SUITE(test_pos_range)

BOOST_AUTO_TEST_CASE(test_pos_range_is_pos_intersect)
{
  // this corresponds to zero-index range [9,19] :
  const pos_range pr(9, 20);

  BOOST_REQUIRE(!pr.is_pos_intersect(8));
  BOOST_REQUIRE(pr.is_pos_intersect(9));
  BOOST_REQUIRE(pr.is_pos_intersect(19));
  BOOST_REQUIRE(!pr.is_pos_intersect(20));
}

BOOST_AUTO_TEST_CASE(test_pos_range_semibound_is_pos_intersect)
{
  // this corresponds to zero-index range [-inf,19] :
  pos_range pr;
  pr.set_end_pos(20);

  BOOST_REQUIRE(pr.is_pos_intersect(8));
  BOOST_REQUIRE(pr.is_pos_intersect(9));
  BOOST_REQUIRE(pr.is_pos_intersect(19));
  BOOST_REQUIRE(!pr.is_pos_intersect(20));
}

BOOST_AUTO_TEST_CASE(test_pos_range_is_range_intersect)
{
  // this corresponds to zero-index range [9,19] :
  const pos_range pr(9, 20);

  // left-side:
  BOOST_REQUIRE(!pr.is_range_intersect(pos_range(0, 9)));
  BOOST_REQUIRE(pr.is_range_intersect(pos_range(0, 10)));

  // right side:
  BOOST_REQUIRE(pr.is_range_intersect(pos_range(19, 30)));
  BOOST_REQUIRE(!pr.is_range_intersect(pos_range(20, 30)));

  // superset:
  BOOST_REQUIRE(pr.is_range_intersect(pos_range(0, 30)));

  // subset:
  BOOST_REQUIRE(pr.is_range_intersect(pos_range(12, 15)));
}

BOOST_AUTO_TEST_CASE(test_pos_range_is_superset_of)
{
  // this corresponds to zero-index range [9,19] :
  const known_pos_range pr(9, 20);

  // non subset tests:
  BOOST_REQUIRE(!pr.is_superset_of(pos_range(7, 8)));
  BOOST_REQUIRE(!pr.is_superset_of(pos_range(8, 10)));
  BOOST_REQUIRE(!pr.is_superset_of(pos_range(8, 20)));
  BOOST_REQUIRE(!pr.is_superset_of(pos_range(9, 21)));
  BOOST_REQUIRE(!pr.is_superset_of(pos_range(8, 21)));
  BOOST_REQUIRE(!pr.is_superset_of(pos_range(25, 30)));

  // subset tests:
  BOOST_REQUIRE(pr.is_superset_of(pos_range(9, 20)));
  BOOST_REQUIRE(pr.is_superset_of(pos_range(11, 20)));
  BOOST_REQUIRE(pr.is_superset_of(pos_range(9, 17)));
  BOOST_REQUIRE(pr.is_superset_of(pos_range(11, 17)));
}

BOOST_AUTO_TEST_SUITE_END()
