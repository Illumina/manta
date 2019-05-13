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
#include "test/testSVLocusScanner.hpp"

#include "SVScorePairRefProcessor.cpp"
#include "SVScorePairRefProcessor.hpp"

BOOST_AUTO_TEST_SUITE(SVScorePairRefProcessor_test_suite)

// Test whether a fragment will support a breakpoint for the ref allele if the following conditions are
// satisfied:
// 1. Start of bam read should overlap with the search range
// 2. Fragment length should be greater than or equal to min fragment length of the sample
// 3. Fragment length should be less than or equal to max fragment length of the sample
// 4. minimum of (breakpoint center pos - fragment start + 1) and (fragment end - breakpoint center pos)
//    should be greater than minimum fragment threshold which is 50
// 5. For RNA, bam read should be properly paired.
// 6. If mate location is less than this read location, fragment support is calculated based on mate's start.
BOOST_AUTO_TEST_CASE(test_processClearedRecord)
{
  const std::vector<bool>         bamFileInfo = {false};  // whether normal bam or tumor bam file
  const bam_header_info           bamHeader(buildTestBamHeader());
  std::unique_ptr<SVLocusScanner> scanner(buildTestSVLocusScanner(bamHeader));
  const PairOptions               options1(false);  // false means DNA options
  SVCandidate                     candidate;
  candidate.insertSeq =
      "GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGT"
      "ATTTTCGTCTGGGGGGTGTGCACGCGATAGCATTGCGAGACGCTGGA";
  // BP1 center pos = 159
  candidate.bp1.interval.range = known_pos_range2(100, 220);
  // BP2 center pos = 309
  candidate.bp2.interval.range = known_pos_range2(250, 370);
  SVEvidence evidence;
  evidence.samples.resize(1);
  SVScorePairRefProcessor processor1(bamFileInfo, scanner.operator*(), options1, candidate, true, evidence);
  // So Search range start = 159 - (125-50) = 84
  // range end = 159 + (125-50) + 1 = 235
  processor1.nextBamIndex(0);
  SVId                       id;
  SVEvidenceWriterSampleData suppFrags;

  std::string querySeq1 = "TCTATCACCCATTTTACCACTCACGGGAGCTCTCC";
  // case-1 is designed here.
  // bam read start = 9. It is not overlapping with search range [84, 235).
  // As a result of this, this fragment is not supporting allele on BP1.
  bam_record bamRecord1;
  buildTestBamRecord(bamRecord1, 0, 9, 0, 100, 35, 15, "35M", querySeq1, 150);
  bamRecord1.set_qname("bamRecord1");
  processor1.processClearedRecord(id, bamRecord1, suppFrags);
  BOOST_REQUIRE(!evidence.getSampleEvidence(0)[bamRecord1.qname()].ref.bp1.isFragmentSupport);

  // case-2 is designed here.
  // Bam read start = 109. It is overlapping with the search range. But
  // its fragment length is 49 which is less than minimum fragment length(50) of the sample
  // As a result of this, this fragment is not supporting allele on BP1.
  bam_record bamRecord2;
  buildTestBamRecord(bamRecord2, 0, 109, 0, 125, 35, 15, "35M", querySeq1, 49);
  bamRecord2.set_qname("bamRecord2");
  processor1.processClearedRecord(id, bamRecord2, suppFrags);
  BOOST_REQUIRE(!evidence.getSampleEvidence(0)[bamRecord2.qname()].ref.bp1.isFragmentSupport);

  // Case-3 is designed here.
  // Bam read start = 109. It is overlapping with the search range. But
  // its fragment length is 130 which is greater than maximum fragment length(125) of the
  // sample. As a result of this, this fragment is not supporting allele on BP1.
  bam_record bamRecord3;
  buildTestBamRecord(bamRecord3, 0, 109, 0, 200, 35, 15, "35M", querySeq1, 130);
  bamRecord3.set_qname("bamRecord3");
  processor1.processClearedRecord(id, bamRecord3, suppFrags);
  BOOST_REQUIRE(!evidence.getSampleEvidence(0)[bamRecord3.qname()].ref.bp1.isFragmentSupport);

  // Case-4 is designed here.
  // Here min(159-109+1, 168-159) = 9 which is less than 50. As a result of this, this fragment
  // is not supporting allele on BP1.
  bam_record bamRecord4;
  buildTestBamRecord(bamRecord4, 0, 109, 0, 125, 35, 15, "35M", querySeq1, 60);
  bamRecord4.set_qname("bamRecord4");
  processor1.processClearedRecord(id, bamRecord4, suppFrags);
  BOOST_REQUIRE(!evidence.getSampleEvidence(0)[bamRecord4.qname()].ref.bp1.isFragmentSupport);

  // Here min(159-109+1, 208-159) = 51 which is greater than 50
  // Fragment start = 109 which is overlapping with the search range[84,235).
  // Fragment size = 100 which is greater than 50 and less than 125.
  // All the above 4 points satisfied, this fragment is supporting allele on BP1.
  bam_record bamRecord5;
  buildTestBamRecord(bamRecord5, 0, 109, 0, 200, 35, 15, "35M", querySeq1, 100);
  bamRecord5.set_qname("bamRecord5");
  processor1.processClearedRecord(id, bamRecord5, suppFrags);
  BOOST_REQUIRE(evidence.getSampleEvidence(0)[bamRecord5.qname()].ref.bp1.isFragmentSupport);

  // The following two test cases are for RNA sample
  const PairOptions       options2(true);
  SVScorePairRefProcessor processor2(bamFileInfo, scanner.operator*(), options2, candidate, true, evidence);
  // For RNA, read should be in proper pair. Following
  // bam read is not in proper pair. As a result of this, this fragment is
  // not supporting allele on BP1.
  bam_record bamRecord6;
  buildTestBamRecord(bamRecord6, 0, 109, 0, 200, 35, 15, "35M", querySeq1, 150);
  bamRecord6.set_qname("bamRecord6");
  processor2.nextBamIndex(0);
  processor2.processClearedRecord(id, bamRecord6, suppFrags);
  BOOST_REQUIRE(!evidence.getSampleEvidence(0)[bamRecord6.qname()].ref.bp1.isFragmentSupport);

  SVScorePairRefProcessor processor3(bamFileInfo, scanner.operator*(), options2, candidate, true, evidence);
  // All the above points are satisfied for RNA sample.
  // Also fragment is in proper pair.
  // So this fragment is supporting allele on BP1.
  bam_record bamRecord7;
  buildTestBamRecord(bamRecord7, 0, 109, 0, 200, 35, 15, "35M", querySeq1, 150);
  bamRecord7.set_qname("bamRecord7");
  bamRecord7.get_data()->core.flag ^= BAM_FLAG::PROPER_PAIR;
  processor3.nextBamIndex(0);
  processor3.processClearedRecord(id, bamRecord7, suppFrags);
  BOOST_REQUIRE(evidence.getSampleEvidence(0)[bamRecord7.qname()].ref.bp1.isFragmentSupport);

  // Case-6 is designed here where mate location
  // is less than this read's start location. Fragment support is
  // calculated based on mate's location and the fragment length.
  // Here min(159-109+1, 208-159) = 51 which is greater than 50
  // Fragment start = 109 which is overlapping with the search range[84,235).
  // Fragment size = 100 which is greater than 50 and less than 125.
  // All the first 4 points satisfied here. As a result of this, the fragment is
  // supporting allele on BP1.
  bam_record bamRecord8;
  buildTestBamRecord(bamRecord8, 0, 180, 0, 109, 84, 15, "84M", "", 100);
  bamRecord8.set_qname("bamRecord8");
  processor1.processClearedRecord(id, bamRecord8, suppFrags);
  BOOST_REQUIRE(evidence.getSampleEvidence(0)[bamRecord8.qname()].ref.bp1.isFragmentSupport);
}

BOOST_AUTO_TEST_SUITE_END()
