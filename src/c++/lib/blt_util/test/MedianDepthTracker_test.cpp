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

#include "MedianDepthTracker.hpp"

BOOST_AUTO_TEST_SUITE(test_MedianDepthTracker)

BOOST_AUTO_TEST_CASE(test_MDT0)
{
  static const double eps(0.00001);

  MedianDepthTracker t;

  BOOST_REQUIRE_CLOSE(t.getMedian(), 0.0, eps);
}

BOOST_AUTO_TEST_CASE(test_MDT1)
{
  static const double eps(0.00001);

  MedianDepthTracker t;

  t.addObs(0);
  t.addObs(2);
  t.addObs(1);
  t.addObs(3);

  BOOST_REQUIRE_CLOSE(t.getMedian(), 2., eps);

  t.addObs(4);

  BOOST_REQUIRE_CLOSE(t.getMedian(), 2.5, eps);
}

BOOST_AUTO_TEST_CASE(test_MDT2)
{
  static const double eps(0.00001);

  MedianDepthTracker t;

  t.addObs(0);
  t.addObs(2);
  t.addObs(1);
  t.addObs(3);

  BOOST_REQUIRE_CLOSE(t.getMedian(), 2., eps);

  t.addObs(2);
  t.addObs(2);
  t.addObs(2);
  t.addObs(2);
  t.addObs(2);

  BOOST_REQUIRE_CLOSE(t.getMedian(), 2., eps);

  t.addObs(1);
  t.addObs(1);
  t.addObs(1);
  t.addObs(1);
  t.addObs(1);
  t.addObs(1);
  t.addObs(1);
  t.addObs(1);
  t.addObs(1);

  BOOST_REQUIRE_CLOSE(t.getMedian(), 1., eps);
}

BOOST_AUTO_TEST_CASE(test_MDT3)
{
  static const double eps(0.00001);

  MedianDepthTracker t;

  t.addObs(1);
  t.addObs(4);
  t.addObs(1);
  t.addObs(4);
  t.addObs(1);
  t.addObs(4);
  t.addObs(1);
  t.addObs(4);
  t.addObs(1);
  t.addObs(4);

  BOOST_REQUIRE_CLOSE(t.getMedian(), 2.5, eps);
}

BOOST_AUTO_TEST_SUITE_END()
