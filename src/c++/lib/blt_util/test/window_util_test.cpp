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

#include "window_util.hpp"

#include <iostream>

BOOST_AUTO_TEST_SUITE(test_window)

BOOST_AUTO_TEST_CASE(test_window)
{
  static const double tol(0.0001);

  window_average wa(3);
  wa.insert(0);
  BOOST_REQUIRE_EQUAL(wa.size(), 1);
  BOOST_REQUIRE_CLOSE(wa.avg(), 0., tol);
  wa.insert(1);
  BOOST_REQUIRE_EQUAL(wa.size(), 2);
  BOOST_REQUIRE_CLOSE(wa.avg(), 0.5, tol);
  wa.insert(2);
  BOOST_REQUIRE_EQUAL(wa.size(), 3);
  BOOST_REQUIRE_CLOSE(wa.avg(), 1., tol);
  wa.insert(3);
  BOOST_REQUIRE_EQUAL(wa.size(), 3);
  BOOST_REQUIRE_CLOSE(wa.avg(), 2., tol);

  wa.insert_null();
  BOOST_REQUIRE_EQUAL(wa.size(), 2);
  BOOST_REQUIRE_CLOSE(wa.avg(), 2.5, tol);
  wa.insert_null();
  BOOST_REQUIRE_EQUAL(wa.size(), 1);
  BOOST_REQUIRE_CLOSE(wa.avg(), 3., tol);
  wa.insert_null();
  BOOST_REQUIRE_EQUAL(wa.size(), 0);

  wa.insert(0);
  BOOST_REQUIRE_EQUAL(wa.size(), 1);
  BOOST_REQUIRE_CLOSE(wa.avg(), 0., tol);
  wa.insert(1);
  BOOST_REQUIRE_EQUAL(wa.size(), 2);
  BOOST_REQUIRE_CLOSE(wa.avg(), 0.5, tol);
  wa.insert(2);
  BOOST_REQUIRE_EQUAL(wa.size(), 3);
  BOOST_REQUIRE_CLOSE(wa.avg(), 1., tol);
  wa.insert(3);
  BOOST_REQUIRE_EQUAL(wa.size(), 3);
  BOOST_REQUIRE_CLOSE(wa.avg(), 2., tol);
}

BOOST_AUTO_TEST_SUITE_END()
