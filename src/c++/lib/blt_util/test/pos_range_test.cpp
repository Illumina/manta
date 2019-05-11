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
