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

#include "RegionSum.hpp"

BOOST_AUTO_TEST_SUITE(test_RegionSum)

BOOST_AUTO_TEST_CASE(RegionSum_test)
{
  RegionSum<unsigned> rs;
  rs.add(known_pos_range2(3, 7), 1u);
  rs.add(known_pos_range2(4, 5), 2u);

  BOOST_REQUIRE_EQUAL(rs.maxVal(), 3u);
}

BOOST_AUTO_TEST_SUITE_END()
