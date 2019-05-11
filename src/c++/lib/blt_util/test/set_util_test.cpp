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

#include "blt_util/set_util.hpp"

BOOST_AUTO_TEST_SUITE(test_set_util)

BOOST_AUTO_TEST_CASE(test_set_subtract)
{
  std::set<int> A;
  A.insert(1);
  A.insert(4);
  A.insert(5);
  A.insert(6);
  A.insert(9);

  std::set<int> B;
  B.insert(2);
  B.insert(5);
  B.insert(6);
  B.insert(7);
  B.insert(10);

  inplaceSetSubtract(A, B);

  BOOST_REQUIRE_EQUAL(B.size(), 3);
  BOOST_REQUIRE_EQUAL(*B.begin(), 2);
}

BOOST_AUTO_TEST_SUITE_END()
