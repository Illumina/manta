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

#include "CircularCounter.hpp"

BOOST_AUTO_TEST_SUITE(test_CircularCounter)

BOOST_AUTO_TEST_CASE(test_CircularCounter1)
{
  CircularCounter cc(3);

  BOOST_CHECK_EQUAL(cc.count(), 0);

  cc.push(true);
  BOOST_CHECK_EQUAL(cc.count(), 1);

  cc.push(false);
  cc.push(false);
  BOOST_CHECK_EQUAL(cc.count(), 1);
  cc.push(false);
  BOOST_CHECK_EQUAL(cc.count(), 0);
}

BOOST_AUTO_TEST_SUITE_END()
