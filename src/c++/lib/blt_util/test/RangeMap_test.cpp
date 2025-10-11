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

#include "blt_util/RangeMap.hpp"

BOOST_AUTO_TEST_SUITE(test_RangeMap)

BOOST_AUTO_TEST_CASE(test_RangeMap)
{
  RangeMap<int, int> rm;

  rm.getRef(2)    = 12;
  rm.getRef(3000) = 13;
  rm.getRef(6000) = 15;

  rm.erase(2);
  rm.getRef(3000) = 3;
  rm.getRef(9000) = 12;
  rm.getRef(6000) = 2;
  rm.erase(9000);

  BOOST_REQUIRE_EQUAL(rm.getConstRef(3000), 3);
  BOOST_REQUIRE_EQUAL(rm.getConstRef(6000), 2);
}

BOOST_AUTO_TEST_CASE(test_RangeMap2)
{
  RangeMap<int, int> rm;

  rm.getRef(10000) = 12;
  rm.getRef(9000)  = 13;
  rm.getRef(7000)  = 15;

  rm.erase(7000);
  BOOST_REQUIRE_EQUAL(rm.getConstRef(9000), 13);
  BOOST_REQUIRE_EQUAL(rm.getConstRef(10000), 12);
  rm.getRef(3000) = 3;
  rm.getRef(9000) = 12;
  rm.getRef(6000) = 2;
  rm.erase(9000);

  BOOST_REQUIRE_EQUAL(rm.getConstRef(3000), 3);
  BOOST_REQUIRE_EQUAL(rm.getConstRef(6000), 2);
  BOOST_REQUIRE_EQUAL(rm.getConstRef(10000), 12);
}

BOOST_AUTO_TEST_CASE(test_rangeMap3)
{
  RangeMap<int, int> rm;

  rm.getRef(0) += 1;
  rm.getRef(20) += 1;
  rm.getRef(21) += 1;
  rm.erase(0);
  rm.erase(20);
  rm.erase(21);
  rm.getRef(0) += 1;
  rm.getRef(20) += 1;
  rm.getRef(21) += 1;

  BOOST_REQUIRE_EQUAL(rm.getConstRef(20), 1);
  BOOST_REQUIRE_EQUAL(rm.getConstRef(21), 1);
}

BOOST_AUTO_TEST_CASE(test_rangeMap_eraseTo)
{
  RangeMap<int, int> rm;

  rm.getRef(5) += 1;
  rm.getRef(20) += 1;
  rm.getRef(21) += 1;

  rm.eraseTo(0);
  BOOST_REQUIRE_EQUAL(rm.getConstRef(5), 1);
  BOOST_REQUIRE_EQUAL(rm.getConstRef(21), 1);

  rm.eraseTo(30);
  BOOST_REQUIRE(rm.empty());

  rm.getRef(5) += 1;
  rm.getRef(20) += 1;
  rm.getRef(21) += 1;

  rm.eraseTo(20);
  BOOST_REQUIRE_EQUAL(rm.getConstRef(21), 1);

  rm.getRef(10030) += 1;
  rm.getRef(10032) += 1;
  rm.getRef(10034) += 1;

  rm.eraseTo(10030);
  BOOST_REQUIRE_EQUAL(rm.getConstRef(10032), 1);
  BOOST_REQUIRE_EQUAL(rm.getConstRef(10034), 1);

  rm.erase(10034);
  BOOST_REQUIRE(!rm.empty());
  rm.erase(10032);
  BOOST_REQUIRE(rm.empty());
}

BOOST_AUTO_TEST_CASE(test_rangeMap_eraseTo2)
{
  RangeMap<int, int> rm;

  rm.getRef(5) += 1;
  rm.getRef(6) += 1;
  rm.getRef(7) += 1;

  rm.eraseTo(5);

  BOOST_REQUIRE_EQUAL(rm.getConstRefDefault(5, 0), 0);
  BOOST_REQUIRE_EQUAL(rm.getConstRef(6), 1);
  BOOST_REQUIRE_EQUAL(rm.getConstRef(7), 1);
}

BOOST_AUTO_TEST_CASE(test_rangeMap_dataSizeBoundary)
{
  RangeMap<int, int> rm(8);

  rm.getRef(0) += 1;
  rm.getRef(8) += 1;
  BOOST_REQUIRE_EQUAL(rm.getConstRef(0), 1);
  BOOST_REQUIRE_EQUAL(rm.getConstRef(8), 1);
}

BOOST_AUTO_TEST_SUITE_END()
