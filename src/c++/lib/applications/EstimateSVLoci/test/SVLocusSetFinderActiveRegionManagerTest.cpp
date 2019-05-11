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

/// \file
/// \author Trevor Ramsay
///

#include "boost/test/unit_test.hpp"

#include "SVLocusSetFinderActiveRegionManager.hpp"

#include "blt_util/depth_buffer.hpp"
#include "svgraph/SVLocusSet.hpp"
#include "test/testAlignmentDataUtil.hpp"
#include "test/testSVLocusUtil.hpp"

static std::shared_ptr<SVLocusSet> getTestSVLocusSet()
{
  bam_header_info bamHeaderInfo(buildTestBamHeader());

  auto set1Ptr(std::make_shared<SVLocusSet>(SVLocusSetOptions(), bamHeaderInfo));
  // build up dummy graph content for set1:
  SVLocus locus1;
  locusAddPair(locus1, 0, 70, 80, 0, 90, 100);

  SVLocus locus2;
  locusAddPair(locus2, 0, 1100, 1110, 0, 1120, 1130);

  SVLocus locus3;
  locusAddPair(locus3, 0, 2100, 2110, 0, 2120, 2130, false, 5);

  set1Ptr->merge(locus1);
  set1Ptr->merge(locus2);
  set1Ptr->merge(locus3);

  return set1Ptr;
}

BOOST_AUTO_TEST_SUITE(SVLocusSetFinderActiveRegionManager_test_suite)

BOOST_AUTO_TEST_CASE(test_SVLocusSetFinderActiveRegionManager)
{
  // Test that SVLocusSetFinder's denoiser removes nodes
  GenomeInterval region(0, 0, 10000);
  auto           set1Ptr(getTestSVLocusSet());

  static const unsigned               denoiseRegionProtectedBorderSize(100);
  SVLocusSetFinderActiveRegionManager regionManager(
      region, set1Ptr, std::make_shared<depth_buffer_compressible>(), denoiseRegionProtectedBorderSize);

  // before
  BOOST_REQUIRE_EQUAL(set1Ptr->totalNodeCount(), 6);

  // Base case, nothing removed
  regionManager.handle_new_pos_value(1);
  BOOST_REQUIRE_EQUAL(set1Ptr->totalNodeCount(), 6);

  // Just before denoieing is triggered
  regionManager.handle_new_pos_value(1098);
  BOOST_REQUIRE_EQUAL(set1Ptr->totalNodeCount(), 6);

  // First case where denoising can trigger on locus1, it is triggered here because of 1000 base minimum
  // denoiseing region size PLUS 100 base border size (using zero-indexed pos values).
  regionManager.handle_new_pos_value(1099);
  BOOST_REQUIRE_EQUAL(set1Ptr->totalNodeCount(), 4);

  // 1000 base chunk size not met, so no cleaning triggered
  regionManager.handle_new_pos_value(2098);
  BOOST_REQUIRE_EQUAL(set1Ptr->totalNodeCount(), 4);

  // 1000 base chunk size met, so locus2 is cleaned
  regionManager.handle_new_pos_value(2099);
  BOOST_REQUIRE_EQUAL(set1Ptr->totalNodeCount(), 2);

  // locus3 not dropped because it is supported by more evidence
  regionManager.handle_new_pos_value(3000);
  BOOST_REQUIRE_EQUAL(set1Ptr->totalNodeCount(), 2);

  // Verification of final state
  regionManager.flush();
  BOOST_REQUIRE_EQUAL(set1Ptr->totalNodeCount(), 2);
}

BOOST_AUTO_TEST_SUITE_END()
