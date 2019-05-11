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

#include "stream_stat.hpp"

BOOST_AUTO_TEST_SUITE(test_stream_stat)

BOOST_AUTO_TEST_CASE(test_stream_stat)
{
  static const double eps(0.00001);

  stream_stat ss;

  ss.add(3.);
  ss.add(4.);
  ss.add(5.);

  BOOST_REQUIRE_EQUAL(ss.size(), 3);
  BOOST_REQUIRE_CLOSE(ss.min(), 3., eps);
  BOOST_REQUIRE_CLOSE(ss.max(), 5., eps);
  BOOST_REQUIRE_CLOSE(ss.mean(), 4., eps);
  BOOST_REQUIRE_CLOSE(ss.sd(), 1., eps);
}

BOOST_AUTO_TEST_SUITE_END()
