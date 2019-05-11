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
