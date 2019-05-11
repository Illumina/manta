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

#include "io_util.hpp"

#include <iomanip>
#include <sstream>

BOOST_AUTO_TEST_SUITE(test_io_util)

BOOST_AUTO_TEST_CASE(test_StreamScoper)
{
  static const double      f(0.123456789);
  static const std::string strf4("0.1235");
  static const std::string strf2("0.12");

  std::ostringstream oss;

  oss << std::fixed << std::setprecision(4);

  oss << f;
  BOOST_REQUIRE_EQUAL(strf4, oss.str());

  oss.str("");

  {
    StreamScoper ss(oss);
    oss << std::fixed << std::setprecision(2);
    oss << f;
    BOOST_REQUIRE_EQUAL(strf2, oss.str());

    oss.str("");
  }

  oss << f;
  BOOST_REQUIRE_EQUAL(strf4, oss.str());
}

BOOST_AUTO_TEST_SUITE_END()
