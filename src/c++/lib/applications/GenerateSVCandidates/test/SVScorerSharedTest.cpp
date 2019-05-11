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
#include "test/testAlignmentDataUtil.hpp"

#include "SVScorerShared.hpp"

BOOST_AUTO_TEST_SUITE(SVScorereSharedTest_test_suite)

// Test the properties of the reads in a fragment
// Following points need to be tested:
// 1. Whether a bam read already scanned or not
// 2. Whether a read is anchored read or tier2 anchored read
// 3. Fragment evidence read size should match with bam read size
BOOST_AUTO_TEST_CASE(test_setReadEvidence)
{
  SVFragmentEvidenceRead svFragmentEvidenceRead1;
  bam_record             bamRecord1;
  buildTestBamRecord(bamRecord1, 0, 200, 0, 300, 35, 50, "35M");
  bamRecord1.set_qname("Read-1");
  // Minimum mapping qulity = 10 which is the threeshold for anchored read
  // Minimum tier2 mapping qulity = 15 which is the threshold for tier2 anchored read
  // Read Mapping quality is 50 which satisfies both the thresholds.
  setReadEvidence(10, 15, bamRecord1, true, svFragmentEvidenceRead1);
  BOOST_REQUIRE(svFragmentEvidenceRead1.isScanned);
  BOOST_REQUIRE(svFragmentEvidenceRead1.mapq == 50);
  BOOST_REQUIRE(svFragmentEvidenceRead1.isAnchored(false));
  BOOST_REQUIRE(svFragmentEvidenceRead1.isAnchored(true));
  BOOST_REQUIRE(svFragmentEvidenceRead1.size == bamRecord1.read_size());

  SVFragmentEvidenceRead svFragmentEvidenceRead2;
  bam_record             bamRecord2;
  buildTestBamRecord(bamRecord2, 0, 200, 0, 300, 35, 12, "35M");
  bamRecord2.set_qname("Read-2");
  // Minimum anchored mapping qulity = 10
  // Minimum tier2 mapping qulity = 15
  // Read Mapping quality = 12 which satisfies only anchored mapping quality threshold.
  setReadEvidence(10, 15, bamRecord2, true, svFragmentEvidenceRead2);
  BOOST_REQUIRE(svFragmentEvidenceRead2.mapq == 12);
  BOOST_REQUIRE(svFragmentEvidenceRead2.isAnchored(false));
  BOOST_REQUIRE(!svFragmentEvidenceRead2.isAnchored(true));

  SVFragmentEvidenceRead svFragmentEvidenceRead3;
  bam_record             bamRecord3;
  buildTestBamRecord(bamRecord3, 0, 200, 0, 300, 35, 12, "35M");
  bamRecord3.set_qname("Read-3");
  // Minimum anchored mapping qulity = 14
  // Minimum tier2 mapping qulity = 10
  // Read Mapping quality = 12 which satisfies only tier2 mapping quality threshold.
  setReadEvidence(14, 10, bamRecord3, true, svFragmentEvidenceRead3);
  BOOST_REQUIRE(svFragmentEvidenceRead3.mapq == 12);
  BOOST_REQUIRE(!svFragmentEvidenceRead3.isAnchored(false));
  BOOST_REQUIRE(svFragmentEvidenceRead3.isAnchored(true));

  SVFragmentEvidenceRead svFragmentEvidenceRead4;
  svFragmentEvidenceRead4.isScanned = true;
  svFragmentEvidenceRead4.mapq      = 70;
  svFragmentEvidenceRead4.setAnchored(true);
  svFragmentEvidenceRead4.setTier2Anchored(false);
  svFragmentEvidenceRead4.size = 150;
  // Read is already scanned so api will not modify any data.
  setReadEvidence(14, 10, bamRecord1, true, svFragmentEvidenceRead4);
  BOOST_REQUIRE(svFragmentEvidenceRead4.isScanned);
  BOOST_REQUIRE(svFragmentEvidenceRead4.mapq == 70);
  BOOST_REQUIRE(svFragmentEvidenceRead4.isAnchored(false));
  BOOST_REQUIRE(!svFragmentEvidenceRead4.isAnchored(true));
  BOOST_REQUIRE(svFragmentEvidenceRead4.size != bamRecord1.read_size());
}
BOOST_AUTO_TEST_SUITE_END()
