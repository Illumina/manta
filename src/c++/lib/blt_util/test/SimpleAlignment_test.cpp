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

#include "blt_util/SimpleAlignment.hpp"

BOOST_AUTO_TEST_SUITE(SimpleAlignment_path)

BOOST_AUTO_TEST_CASE(test_SimpleAlignment)
{
  using namespace ALIGNPATH;

  SimpleAlignment al;

  al.pos = 100;
  cigar_to_apath("100M", al.path);

  const known_pos_range2 testRange(matchifyEdgeSoftClipRefRange(al));

  BOOST_REQUIRE_EQUAL(testRange.begin_pos(), 100);
  BOOST_REQUIRE_EQUAL(testRange.end_pos(), 200);
}

BOOST_AUTO_TEST_CASE(test_SimpleAlignment2)
{
  using namespace ALIGNPATH;

  SimpleAlignment al;

  al.pos = 100;
  cigar_to_apath("10S50M10D40M10S", al.path);

  const known_pos_range2 testRange(matchifyEdgeSoftClipRefRange(al));

  BOOST_REQUIRE_EQUAL(testRange.begin_pos(), 90);
  BOOST_REQUIRE_EQUAL(testRange.end_pos(), 210);
}

BOOST_AUTO_TEST_SUITE_END()
