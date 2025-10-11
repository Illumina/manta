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
