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

#include "blt_util/set_util.hpp"

BOOST_AUTO_TEST_SUITE(test_set_util)

BOOST_AUTO_TEST_CASE(test_set_subtract)
{
  std::set<int> A;
  A.insert(1);
  A.insert(4);
  A.insert(5);
  A.insert(6);
  A.insert(9);

  std::set<int> B;
  B.insert(2);
  B.insert(5);
  B.insert(6);
  B.insert(7);
  B.insert(10);

  inplaceSetSubtract(A, B);

  BOOST_REQUIRE_EQUAL(B.size(), 3);
  BOOST_REQUIRE_EQUAL(*B.begin(), 2);
}

BOOST_AUTO_TEST_SUITE_END()
