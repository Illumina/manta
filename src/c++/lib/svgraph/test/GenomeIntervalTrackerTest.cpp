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
