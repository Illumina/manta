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

#include "string_util.hpp"

#include <cstring>

BOOST_AUTO_TEST_SUITE(string_util)

static const char* test_string("234562342");

template <typename T>
void test_split_string_bytype(T test)
{
  std::vector<std::string> result;

  split_string(test, '2', result);
  BOOST_REQUIRE_EQUAL(static_cast<int>(result.size()), 4);
  BOOST_REQUIRE_EQUAL(result[0], "");
  BOOST_REQUIRE_EQUAL(result[1], "3456");
  BOOST_REQUIRE_EQUAL(result[3], "");

  split_string(test, '2', result, true);
  BOOST_REQUIRE_EQUAL(static_cast<int>(result.size()), 2);
  BOOST_REQUIRE_EQUAL(result[0], "3456");

  split_string(test, 'X', result);
  BOOST_REQUIRE_EQUAL(static_cast<int>(result.size()), 1);
  BOOST_REQUIRE_EQUAL(result[0], test);

  split_string("", 'X', result);
  BOOST_REQUIRE_EQUAL(static_cast<int>(result.size()), 1);
  BOOST_REQUIRE_EQUAL(result[0], "");
}

BOOST_AUTO_TEST_CASE(test_split_string_cstr)
{
  test_split_string_bytype(test_string);
}

BOOST_AUTO_TEST_CASE(test_split_string)
{
  const std::string test(test_string);
  test_split_string_bytype(test);
}

BOOST_AUTO_TEST_CASE(test_destructive_split_string)
{
  std::unique_ptr<char[]> tmp(new char[strlen(test_string) + 1]);
  std::strcpy(tmp.get(), test_string);

  std::vector<const char*> result;
  destructive_split_string(tmp.get(), '2', result);
  BOOST_REQUIRE_EQUAL(static_cast<int>(result.size()), 4);
  BOOST_REQUIRE_EQUAL(std::strcmp(result[0], ""), 0);
  BOOST_REQUIRE_EQUAL(std::strcmp(result[1], "3456"), 0);
  BOOST_REQUIRE_EQUAL(std::strcmp(result[3], ""), 0);
}

BOOST_AUTO_TEST_CASE(test_split_match)
{
  const std::string test(test_string);

  BOOST_REQUIRE(split_match(test, '2', "34"));
  BOOST_REQUIRE(!split_match(test, '2', "XX"));
}

BOOST_AUTO_TEST_SUITE_END()
