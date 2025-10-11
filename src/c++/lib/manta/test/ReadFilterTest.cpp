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

///
/// \file
/// \author Trevor Ramsay
///

#include "boost/test/unit_test.hpp"

#include "ReadFilter.hpp"

#include "blt_util/align_path.hpp"
#include "htsapi/align_path_bam_util.hpp"
#include "test/testAlignmentDataUtil.hpp"

BOOST_AUTO_TEST_SUITE(test_ReadFilter)

// Test all of the paths possible for core read filtering.
BOOST_AUTO_TEST_CASE(test_isReadFilteredCore)
{
  BOOST_TEST_MESSAGE("SDS MANTA-663");

  ALIGNPATH::path_t inputPath;
  cigar_to_apath("100M2000D100M", inputPath);

  bam_record bamRead;
  bam1_t*    bamDataPtr(bamRead.get_data());
  edit_bam_cigar(inputPath, *bamDataPtr);

  bam_record primaryFilterRead(bamRead);
  primaryFilterRead.toggle_is_filtered();

  bam_record unmappedRead(bamRead);
  unmappedRead.toggle_is_unmapped();

  bam_record duplicateRead(bamRead);
  duplicateRead.toggle_is_duplicate();

  bam_record supplementNoSASplitRead(bamRead);
  supplementNoSASplitRead.toggle_is_supplementary();
  BOOST_REQUIRE(!supplementNoSASplitRead.isSASplit());

  bam_record supplementSASplitRead(bamRead);
  supplementSASplitRead.toggle_is_supplementary();
  addSupplementaryAlignmentEvidence(supplementSASplitRead);
  BOOST_REQUIRE(supplementSASplitRead.isSASplit());

  bam_record secondarySASplitRead(bamRead);
  secondarySASplitRead.toggle_is_secondary();
  addSupplementaryAlignmentEvidence(secondarySASplitRead);
  BOOST_REQUIRE(secondarySASplitRead.isSASplit());

  bam_record secondaryNoSASplitRead(bamRead);
  secondaryNoSASplitRead.toggle_is_secondary();
  BOOST_REQUIRE(!secondaryNoSASplitRead.isSASplit());

  // Test multiple read states for isReadFilteredCore.
  BOOST_REQUIRE(!isReadFilteredCore(bamRead));
  BOOST_REQUIRE(!isReadFilteredCore(unmappedRead));
  BOOST_REQUIRE(isReadFilteredCore(primaryFilterRead));
  BOOST_REQUIRE(!isReadFilteredCore(secondarySASplitRead));
  BOOST_REQUIRE(isReadFilteredCore(secondaryNoSASplitRead));
  BOOST_REQUIRE(!isReadFilteredCore(supplementSASplitRead));
  BOOST_REQUIRE(isReadFilteredCore(supplementNoSASplitRead));
}

// Test the mapped core read filtering.
BOOST_AUTO_TEST_CASE(test_isMappedReadFilteredCore)
{
  BOOST_TEST_MESSAGE("SDS MANTA-663");

  bam_record bamRead;
  buildTestBamRecord(bamRead);

  BOOST_REQUIRE(!isReadUnmappedOrFilteredCore(bamRead));
  bamRead.toggle_is_unmapped();
  BOOST_REQUIRE(isReadUnmappedOrFilteredCore(bamRead));
}

BOOST_AUTO_TEST_SUITE_END()
