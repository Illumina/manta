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

#include "blt_util/seq_util.hpp"

BOOST_AUTO_TEST_SUITE(test_seq_util)

BOOST_AUTO_TEST_CASE(test_seq_repeat)
{
  std::string ru;
  unsigned    count;

  static const std::string test0 = "";
  get_seq_repeat_unit(test0, ru, count);
  BOOST_REQUIRE_EQUAL(ru, "");
  BOOST_REQUIRE_EQUAL(count, 1u);

  static const std::string test1 = "AAAA";
  get_seq_repeat_unit(test1, ru, count);
  BOOST_REQUIRE_EQUAL(ru, "A");
  BOOST_REQUIRE_EQUAL(count, 4u);

  static const std::string test2 = "ACAC";
  get_seq_repeat_unit(test2, ru, count);
  BOOST_REQUIRE_EQUAL(ru, "AC");
  BOOST_REQUIRE_EQUAL(count, 2u);

  static const std::string test3 = "TACAC";
  get_seq_repeat_unit(test3, ru, count);
  BOOST_REQUIRE_EQUAL(ru, "TACAC");
  BOOST_REQUIRE_EQUAL(count, 1u);

  get_vcf_seq_repeat_unit(test3, ru, count);
  BOOST_REQUIRE_EQUAL(ru, "AC");
  BOOST_REQUIRE_EQUAL(count, 2u);
}

BOOST_AUTO_TEST_SUITE_END()
