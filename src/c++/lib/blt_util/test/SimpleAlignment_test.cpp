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
