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

#include "CircularCounter.hpp"

BOOST_AUTO_TEST_SUITE(test_CircularCounter)

BOOST_AUTO_TEST_CASE(test_CircularCounter1)
{
  CircularCounter cc(3);

  BOOST_CHECK_EQUAL(cc.count(), 0);

  cc.push(true);
  BOOST_CHECK_EQUAL(cc.count(), 1);

  cc.push(false);
  cc.push(false);
  BOOST_CHECK_EQUAL(cc.count(), 1);
  cc.push(false);
  BOOST_CHECK_EQUAL(cc.count(), 0);
}

BOOST_AUTO_TEST_SUITE_END()
