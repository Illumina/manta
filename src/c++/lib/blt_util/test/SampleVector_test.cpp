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

#include "SampleVector.hpp"

#include <random>

BOOST_AUTO_TEST_SUITE(test_SampleVector)

BOOST_AUTO_TEST_CASE(test_SampleVector1)
{
  std::mt19937                    rngEngine(0);
  SampleVector<int, std::mt19937> sv(2, rngEngine);

  for (unsigned i(0); i < 100; ++i) {
    sv.push(i);
  }

  // can't make this portable:
  // BOOST_CHECK_EQUAL(sv.data()[0],34);
  // BOOST_CHECK_EQUAL(sv.data()[1],70);
}

BOOST_AUTO_TEST_SUITE_END()
