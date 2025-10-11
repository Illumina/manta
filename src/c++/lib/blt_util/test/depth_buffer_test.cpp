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

#include "depth_buffer.hpp"

BOOST_AUTO_TEST_SUITE(test_depth_buffer)

/// return buffer loaded with simple test pattern
///
/// at depth 10Y, depth is Y
///
static depth_buffer get_db_test_pattern()
{
  depth_buffer db;

  // load a depth pattern in:
  for (unsigned i(100); i < 110; ++i) {
    for (unsigned j(i + 1); j < 110; ++j) {
      db.inc(j);
    }
  }
  return db;
}

/// return buffer loaded with simple test pattern
///
/// at pos Y, depth is Y-100 before compression
///
static depth_buffer_compressible get_db_compressible_test_pattern(const unsigned compressionLevel)
{
  depth_buffer_compressible db(compressionLevel);

  // load a depth pattern in:
  for (unsigned i(101); i < 200; ++i) {
    db.inc(i, 200 - i);
  }
  return db;
}

BOOST_AUTO_TEST_CASE(test_depth_buffer_val)
{
  depth_buffer db(get_db_test_pattern());
  BOOST_REQUIRE_EQUAL(static_cast<int>(db.val(109)), 9);
}

BOOST_AUTO_TEST_CASE(test_depth_buffer_clear)
{
  depth_buffer db(get_db_test_pattern());
  db.clear_pos(109);
  BOOST_REQUIRE_EQUAL(static_cast<int>(db.val(109)), 0);
}

BOOST_AUTO_TEST_CASE(test_depth_buffer_range)
{
  depth_buffer db(get_db_test_pattern());
  BOOST_REQUIRE(!db.is_range_ge_than(0, 107, 8));
  BOOST_REQUIRE(db.is_range_ge_than(0, 108, 8));
}

BOOST_AUTO_TEST_CASE(test_depth_buffer_compressible_val)
{
  depth_buffer_compressible db(get_db_compressible_test_pattern(8));
  BOOST_REQUIRE_EQUAL(static_cast<int>(db.val(107)), 8);
  BOOST_REQUIRE_EQUAL(static_cast<int>(db.val(110)), 8);
  BOOST_REQUIRE_EQUAL(static_cast<int>(db.val(150)), 48);
  BOOST_REQUIRE_EQUAL(static_cast<int>(db.val(199)), 96);
}

BOOST_AUTO_TEST_CASE(test_depth_buffer_compressible_clear)
{
  depth_buffer_compressible db(get_db_compressible_test_pattern(8));
  for (unsigned i(100); i < 120; ++i) {
    db.clear_pos(i);
  }
  BOOST_REQUIRE_EQUAL(static_cast<int>(db.val(110)), 0);
  BOOST_REQUIRE_EQUAL(static_cast<int>(db.val(150)), 48);
  BOOST_REQUIRE_EQUAL(static_cast<int>(db.val(199)), 96);
  BOOST_REQUIRE_EQUAL(static_cast<int>(db.val(109)), 0);
}

BOOST_AUTO_TEST_SUITE_END()
