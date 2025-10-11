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

#include "id_map.hpp"

#include <string>

BOOST_AUTO_TEST_SUITE(test_id_map)

BOOST_AUTO_TEST_CASE(test_id_set)
{
  id_set<std::string> iset;

  iset.insert_key("brown");
  iset.insert_key("fox");

  BOOST_REQUIRE_EQUAL(iset.test_key("brown"), true);
  BOOST_REQUIRE_EQUAL(iset.test_key("x"), false);

  const unsigned expect_id(1);

  BOOST_REQUIRE_EQUAL(iset.get_id("fox"), expect_id);
  BOOST_REQUIRE_EQUAL(iset.get_key(expect_id), std::string("fox"));
}

BOOST_AUTO_TEST_CASE(test_id_map)
{
  id_map<std::string, std::string> imap;

  imap.insert("brown", "123");
  imap.insert("fox", "456");

  BOOST_REQUIRE_EQUAL(imap.test_key("brown"), true);
  BOOST_REQUIRE_EQUAL(imap.test_key("x"), false);

  // test key replacement:
  imap.insert("brown", "789");
  BOOST_REQUIRE_EQUAL(imap.get_value(imap.get_id("brown")), std::string("789"));

  const unsigned expect_id(1);

  BOOST_REQUIRE_EQUAL(imap.get_id("fox"), expect_id);
  BOOST_REQUIRE_EQUAL(imap.get_key(expect_id), std::string("fox"));
  BOOST_REQUIRE_EQUAL(imap.get_value(expect_id), std::string("456"));
}

BOOST_AUTO_TEST_SUITE_END()
