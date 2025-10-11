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
