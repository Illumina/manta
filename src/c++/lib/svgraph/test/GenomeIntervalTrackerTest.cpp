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

#include "svgraph/GenomeIntervalTracker.hpp"

BOOST_AUTO_TEST_SUITE(GenomeIntervalTracker_tests)

BOOST_AUTO_TEST_CASE(GenomIntervalTracker_subset)
{
  GenomeIntervalTracker git;

  git.addInterval(GenomeInterval(1, 100, 200));
  git.addInterval(GenomeInterval(1, 150, 250));
  git.addInterval(GenomeInterval(3, 100, 200));

  BOOST_REQUIRE(git.isSubsetOfRegion(GenomeInterval(1, 170, 220)));
  BOOST_REQUIRE(!git.isSubsetOfRegion(GenomeInterval(1, 220, 280)));
  BOOST_REQUIRE(!git.isSubsetOfRegion(GenomeInterval(1, 300, 350)));
  BOOST_REQUIRE(!git.isSubsetOfRegion(GenomeInterval(2, 150, 200)));
}

BOOST_AUTO_TEST_CASE(GenomIntervalTracker_clear)
{
  GenomeIntervalTracker git;

  git.addInterval(GenomeInterval(1, 100, 200));
  git.addInterval(GenomeInterval(1, 150, 250));
  git.addInterval(GenomeInterval(3, 100, 200));

  git.clear();

  BOOST_REQUIRE(!git.isSubsetOfRegion(GenomeInterval(1, 150, 200)));
}

BOOST_AUTO_TEST_CASE(GenomIntervalTracker_merge)
{
  GenomeIntervalTracker git1, git2;

  git1.addInterval(GenomeInterval(1, 100, 200));
  git2.addInterval(GenomeInterval(1, 150, 250));
  git2.addInterval(GenomeInterval(3, 100, 200));

  git1.merge(git2);

  BOOST_REQUIRE(git1.isSubsetOfRegion(GenomeInterval(1, 170, 220)));
  BOOST_REQUIRE(git1.isSubsetOfRegion(GenomeInterval(3, 170, 180)));
  BOOST_REQUIRE(!git1.isSubsetOfRegion(GenomeInterval(1, 220, 280)));
  BOOST_REQUIRE(!git1.isSubsetOfRegion(GenomeInterval(1, 300, 350)));
  BOOST_REQUIRE(!git1.isSubsetOfRegion(GenomeInterval(2, 150, 200)));
}

BOOST_AUTO_TEST_SUITE_END()
