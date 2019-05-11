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

// hack to call private methods of SVLocusSet:
//#pragma clang diagnostic ignored "-Wkeyword-macro"
//#define private public

#include "svgraph/SVLocusSet.hpp"
#include "test/testSVLocusUtil.hpp"

static unsigned testOverlap(
    SVLocusSet& locusSet, const int32_t tid, const int32_t beginPos, const int32_t endPos)
{
  std::set<SVLocusSet::NodeAddressType> intersect;
  locusSet.getRegionIntersect(GenomeInterval(tid, beginPos, endPos), intersect);
  return intersect.size();
}

BOOST_AUTO_TEST_SUITE(SVLocusSetPrivate_test_suite)

BOOST_AUTO_TEST_CASE(test_SVLocusIntersect)
{
  // construct a simple two-node locus
  SVLocus locus1;
  locusAddPair(locus1, 1, 10, 20, 2, 30, 40);

  SVLocusSet set1;
  set1.merge(locus1);
  set1.checkState(true, true);

  // test for various intersections:

  // non-overlap test:
  BOOST_REQUIRE_EQUAL(testOverlap(set1, 1, 1, 2), 0u);

  // left-edge overlap:
  BOOST_REQUIRE_EQUAL(testOverlap(set1, 1, 9, 11), 1u);

  // right-edge overlap:
  BOOST_REQUIRE_EQUAL(testOverlap(set1, 1, 19, 21), 1u);

  // non-overlap:
  BOOST_REQUIRE_EQUAL(testOverlap(set1, 1, 29, 31), 0u);

  // non-overlap (diff tid):
  BOOST_REQUIRE_EQUAL(testOverlap(set1, 2, 9, 11), 0u);
}

BOOST_AUTO_TEST_CASE(test_SVLocusCombine)
{
  // test reassigning the locus numbers of non-overlapping loci in a set:

  // construct a simple two-node locus
  SVLocus locus1;
  locusAddPair(locus1, 1, 10, 20, 2, 30, 40);

  SVLocus locus2;
  locusAddPair(locus2, 3, 10, 20, 4, 30, 40);

  SVLocus locus3;
  locusAddPair(locus3, 5, 10, 20, 6, 30, 40);

  SVLocusSetOptions sopt;
  sopt.minMergeEdgeObservations = 1;
  SVLocusSet set1(sopt);
  set1.merge(locus1);
  set1.merge(locus2);
  set1.merge(locus3);

  set1.checkState(true, true);

  const SVLocusSet& cset1(set1);

  BOOST_REQUIRE_EQUAL(cset1.getLocus(0).size(), 2u);
  BOOST_REQUIRE_EQUAL(cset1.getLocus(1).size(), 2u);
  BOOST_REQUIRE_EQUAL(cset1.getLocus(2).size(), 2u);
#if 0
    set1.combineLoci(0,0);
    set1.combineLoci(2,0);

    BOOST_REQUIRE_EQUAL(cset1.getLocus(0).size(),4u);
    BOOST_REQUIRE_EQUAL(cset1.getLocus(1).size(),2u);
    BOOST_REQUIRE_EQUAL(cset1.getLocus(2).size(),0u);
#endif
}

BOOST_AUTO_TEST_SUITE_END()
