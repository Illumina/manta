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
