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

#include "blt_util/known_pos_range2.hpp"

BOOST_AUTO_TEST_SUITE(test_known_pos_range2)

BOOST_AUTO_TEST_CASE(test_known_pos_range2_is_pos_intersect)
{
  // this corresponds to zero-index range [9,19] :
  const known_pos_range2 pr(9, 20);

  BOOST_REQUIRE(!pr.is_pos_intersect(8));
  BOOST_REQUIRE(pr.is_pos_intersect(9));
  BOOST_REQUIRE(pr.is_pos_intersect(19));
  BOOST_REQUIRE(!pr.is_pos_intersect(20));
}

BOOST_AUTO_TEST_CASE(test_pos_range_semibound_is_pos_intersect)
{
  // this corresponds to zero-index range [-inf,19] :
  known_pos_range2 pr;
  pr.set_end_pos(20);

  BOOST_REQUIRE(pr.is_pos_intersect(8));
  BOOST_REQUIRE(pr.is_pos_intersect(9));
  BOOST_REQUIRE(pr.is_pos_intersect(19));
  BOOST_REQUIRE(!pr.is_pos_intersect(20));
}

BOOST_AUTO_TEST_CASE(test_known_pos_range2_is_range_intersect)
{
  // this corresponds to zero-index range [9,19] :
  const known_pos_range2 pr(9, 20);

  // left-side:
  BOOST_REQUIRE(!pr.is_range_intersect(known_pos_range2(0, 9)));
  BOOST_REQUIRE(pr.is_range_intersect(known_pos_range2(0, 10)));

  // right side:
  BOOST_REQUIRE(pr.is_range_intersect(known_pos_range2(19, 30)));
  BOOST_REQUIRE(!pr.is_range_intersect(known_pos_range2(20, 30)));

  // superset:
  BOOST_REQUIRE(pr.is_range_intersect(known_pos_range2(0, 30)));

  // subset:
  BOOST_REQUIRE(pr.is_range_intersect(known_pos_range2(12, 15)));
}

BOOST_AUTO_TEST_CASE(test_known_pos_range2_intersect_window)
{
  // this corresponds to zero-index range [9,19] :
  const known_pos_range2 pr1(9, 20);
  const known_pos_range2 pr2(30, 40);

  BOOST_REQUIRE(!is_intersect_window(pr1, pr2));
  BOOST_REQUIRE(!is_intersect_window(pr2, pr1));

  BOOST_REQUIRE(!is_intersect_window(pr1, pr2, 10));
  BOOST_REQUIRE(!is_intersect_window(pr2, pr1, 10));

  BOOST_REQUIRE(is_intersect_window(pr1, pr2, 11));
  BOOST_REQUIRE(is_intersect_window(pr2, pr1, 11));
}

BOOST_AUTO_TEST_SUITE_END()
