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

#include "boost/test/floating_point_comparison.hpp"
#include "boost/test/unit_test.hpp"

#include "blt_util/SizeDistribution.cpp"

BOOST_AUTO_TEST_SUITE(test_SizeDistribution)

BOOST_AUTO_TEST_CASE(test_EmptySizeDistribution)
{
  SizeDistribution sd;

  BOOST_REQUIRE_EQUAL(sd.cdf(2), 0);
  BOOST_REQUIRE_EQUAL(sd.quantile(0.2), 0);
}

BOOST_AUTO_TEST_CASE(test_SizeDistribution1)
{
  SizeDistribution sd;

  sd.addObservation(1);
  sd.addObservation(2);
  sd.addObservation(3);
  sd.addObservation(4);

  BOOST_REQUIRE_EQUAL(sd.cdf(0), 0.);
  BOOST_REQUIRE_EQUAL(sd.cdf(1), 0.25);
  BOOST_REQUIRE_EQUAL(sd.cdf(2), 0.5);
  BOOST_REQUIRE_EQUAL(sd.cdf(3), 0.75);
  BOOST_REQUIRE_EQUAL(sd.cdf(4), 1);
  BOOST_REQUIRE_EQUAL(sd.quantile(0.0), 1);
  BOOST_REQUIRE_EQUAL(sd.quantile(0.25), 1);
  BOOST_REQUIRE_EQUAL(sd.quantile(0.5), 2);
  BOOST_REQUIRE_EQUAL(sd.quantile(0.74), 3);
  BOOST_REQUIRE_EQUAL(sd.quantile(0.75), 3);
  BOOST_REQUIRE_EQUAL(sd.quantile(0.76), 4);
  BOOST_REQUIRE_EQUAL(sd.quantile(1.0), 4);
}

BOOST_AUTO_TEST_CASE(test_SizeDistributionFilter)
{
  SizeDistribution sd;

  sd.addObservation(1);
  sd.addObservation(2);
  sd.addObservation(3);
  sd.addObservation(4);

  sd.filterObservationsOverQuantile(0.5);

  BOOST_REQUIRE_EQUAL(sd.totalObservations(), 2u);

  BOOST_REQUIRE_EQUAL(sd.cdf(0), 0.);
  BOOST_REQUIRE_EQUAL(sd.cdf(1), 0.5);
  BOOST_REQUIRE_EQUAL(sd.cdf(2), 1);
  BOOST_REQUIRE_EQUAL(sd.quantile(0.0), 1);
  BOOST_REQUIRE_EQUAL(sd.quantile(0.25), 1);
  BOOST_REQUIRE_EQUAL(sd.quantile(0.5), 1);
  BOOST_REQUIRE_EQUAL(sd.quantile(0.75), 2);
  BOOST_REQUIRE_EQUAL(sd.quantile(1.0), 2);
}

BOOST_AUTO_TEST_CASE(test_SizeDistributionPdf)
{
  SizeDistribution sd;

  sd.addObservation(1);
  sd.addObservation(2);
  sd.addObservation(3);
  sd.addObservation(4);
  sd.addObservation(5);
  sd.addObservation(6);
  sd.addObservation(7);
  sd.addObservation(8);
  sd.addObservation(9);
  sd.addObservation(10);

  static const float tol(0.0001);
  BOOST_REQUIRE_CLOSE(sd.pdf(1), 0.1f, tol);
  BOOST_REQUIRE_CLOSE(sd.pdf(5), 0.1f, tol);
  BOOST_REQUIRE_CLOSE(sd.pdf(10), 0.1f, tol);

  static const float expect2(0.1 * 5. / 6.);
  BOOST_REQUIRE_CLOSE(sd.pdf(0), expect2, tol);
  BOOST_REQUIRE_CLOSE(sd.pdf(11), expect2, tol);
}

BOOST_AUTO_TEST_SUITE_END()
