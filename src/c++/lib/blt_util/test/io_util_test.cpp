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
