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
