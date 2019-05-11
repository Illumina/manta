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
