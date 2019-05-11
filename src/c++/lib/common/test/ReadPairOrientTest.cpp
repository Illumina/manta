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

#include "ReadPairOrient.hpp"

BOOST_AUTO_TEST_SUITE(test_ReadPairOrient)

BOOST_AUTO_TEST_CASE(test_PairTypes)
{
  // innies:
  {
    const pos_t readAPos(10);
    const bool  readAFwd(true);
    const pos_t readBPos(20);
    const bool  readBFwd(false);

    // map A->1 B->2
    PAIR_ORIENT::index_t res = PAIR_ORIENT::get_index(readAPos, readAFwd, readBPos, readBFwd);
    BOOST_REQUIRE_EQUAL(res, PAIR_ORIENT::Rp);

    // map A->2 B->1
    PAIR_ORIENT::index_t res2 = PAIR_ORIENT::get_index(readBPos, readBFwd, readAPos, readAFwd);
    BOOST_REQUIRE_EQUAL(res2, PAIR_ORIENT::Rp);
  }

  // outties
  {
    const pos_t readAPos(30);
    const bool  readAFwd(true);
    const pos_t readBPos(20);
    const bool  readBFwd(false);

    // map A->1 B->2
    PAIR_ORIENT::index_t res = PAIR_ORIENT::get_index(readAPos, readAFwd, readBPos, readBFwd);
    BOOST_REQUIRE_EQUAL(res, PAIR_ORIENT::Rm);

    // map A->2 B->1
    PAIR_ORIENT::index_t res2 = PAIR_ORIENT::get_index(readBPos, readBFwd, readAPos, readAFwd);
    BOOST_REQUIRE_EQUAL(res2, PAIR_ORIENT::Rm);
  }

  // short fragments should resolve to innies:
  {
    const pos_t readAPos(10);
    const bool  readAFwd(true);
    const pos_t readBPos(10);
    const bool  readBFwd(false);

    // map A->1 B->2
    PAIR_ORIENT::index_t res = PAIR_ORIENT::get_index(readAPos, readAFwd, readBPos, readBFwd);
    BOOST_REQUIRE_EQUAL(res, PAIR_ORIENT::Rp);

    // map A->2 B->1
    PAIR_ORIENT::index_t res2 = PAIR_ORIENT::get_index(readBPos, readBFwd, readAPos, readAFwd);
    BOOST_REQUIRE_EQUAL(res2, PAIR_ORIENT::Rp);
  }
}

BOOST_AUTO_TEST_SUITE_END()
