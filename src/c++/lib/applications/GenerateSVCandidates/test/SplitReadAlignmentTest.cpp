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
/// \author Atanu Pal
///

#include <iostream>

#include "boost/test/unit_test.hpp"

#include "options/CallOptionsShared.hpp"
#include "test/testAlignmentDataUtil.hpp"

#include "SplitReadAlignment.cpp"
#include "SplitReadAlignment.hpp"

BOOST_AUTO_TEST_SUITE(SplitReadAligment_test_suite)

// Check evidence based on the alignment information
// Consider the schematic diagram below:
//                                |<-  BP  ->|
// Alignment->    --------------------------------------------->
//                |    left size  | hom size |   right size |
// Following Points need to tested:
// 1. Left size should be greater than equal to minFlank size
// 2. Right size should be greater than equal to minFlank size
// 3. Fraction of mismatches in left side less than 0.25
// 4. Fraction of mismatches in right side less than 0.25
// 5. Fraction of alignment score should be greater than 0.9 where
// alignment score = Total Size - (Left side mismatches + Hom side mismatches + Right side mismatches)
BOOST_AUTO_TEST_CASE(test_ISEvidenceCheck)
{
  SRAlignmentInfo alignmentInfo;

  alignmentInfo.leftMismatches  = 1;
  alignmentInfo.rightMismatches = 1;
  alignmentInfo.leftSize        = 15;
  alignmentInfo.rightSize       = 15;
  alignmentInfo.alignScore      = 28;
  // Here all values are correct
  // min flank size = 5
  // Fraction of mismatches in left side = 0.067(1/15)
  // Fraction of mismatches in right side = 0.067(1/15)
  // Fraction of alignment score = 0.9333 (28/30)
  BOOST_REQUIRE(isEvidenceCheck(alignmentInfo, 5));

  alignmentInfo.leftSize = 10;
  // min flank size = 11. API should return false as left size
  // is less than 11.
  BOOST_REQUIRE(!isEvidenceCheck(alignmentInfo, 11));

  alignmentInfo.rightSize = 8;
  // min flank size = 9. Althouh left size = 10 (>minFlankSize), then also API should return false
  // as right size is less than 9.
  BOOST_REQUIRE(!isEvidenceCheck(alignmentInfo, 9));

  alignmentInfo.leftMismatches = 3;
  // fraction of mismatches in left side = 0.3(3/10) which is greater than 0.25.
  BOOST_REQUIRE(!isEvidenceCheck(alignmentInfo, 5));

  alignmentInfo.leftMismatches  = 1;
  alignmentInfo.rightMismatches = 3;
  // fraction of mismatches in left side = 0.1 (1/10) which is less than 0.25, but
  // fraction of mismatches in right side = 0.375 (3/8) which is greater than 0.25.
  BOOST_REQUIRE(!isEvidenceCheck(alignmentInfo, 5));

  alignmentInfo.leftMismatches  = 1;
  alignmentInfo.rightMismatches = 1;
  alignmentInfo.alignScore      = 15;
  // Here fraction of mismatches satisfies the criteria but
  // fraction of alignment score (0.833 = 15/18) which is less
  // than 0.9.
  BOOST_REQUIRE(!isEvidenceCheck(alignmentInfo, 5));
}

// Test the evidence value of an alignment which is
// 2*min(left size, right size) / (left size + right size). Left size and right
// size are described in test_ISEvidenceCheck.
// If an alignment satisfies test_ISEvidenceCheck either for minFlankSize=16 or
// minFlankSizeTier2=8, then the non-zero evidence value will be set.
BOOST_AUTO_TEST_CASE(test_setEvidence)
{
  SRAlignmentInfo alignmentInfo;

  alignmentInfo.leftSize  = 10;
  alignmentInfo.rightSize = 7;
  // alignmentInfo is not satisfied for test_ISEvidenceCheck as minFlankSize=16
  // So evidence is 0 for this.
  setEvidence(alignmentInfo);
  BOOST_REQUIRE_EQUAL(alignmentInfo.evidence, 0);

  alignmentInfo.leftSize   = 20;
  alignmentInfo.rightSize  = 20;
  alignmentInfo.alignScore = 40;
  // alignmentInfo is  satisfied for test_ISEvidenceCheck as minFlankSize=16
  // So evidence is 1 for this.
  setEvidence(alignmentInfo);
  BOOST_REQUIRE_EQUAL(alignmentInfo.evidence, 1);

  alignmentInfo.leftSize   = 10;
  alignmentInfo.rightSize  = 10;
  alignmentInfo.alignScore = 20;
  // alignmentInfo is  not satisfied for test_ISEvidenceCheck for minFlankSize=16,
  // but it is satisfied for minFlankSizeTier2=8.
  // So evidence is 1 for this.
  setEvidence(alignmentInfo);
  BOOST_REQUIRE_EQUAL(alignmentInfo.evidence, 1);
}

// Test the alignment score of an alignment
// Alignment score = total size - (left mismatch + hom mismatch + right mismatch)
// left mismatch, hom mismatch and right mismatch are mentioned in test_ISEvidenceCheck.
BOOST_AUTO_TEST_CASE(test_calculateAlignScore)
{
  std::string targetSeq =
      "GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGT"
      "ATTTTCGTCTGGGGGGTGTGCACGCGATAGCATTGCGAGACGCTGGA";

  // Alignment of following query sequence starts at position 9 (0-based)
  // with 10 matches in the left side followed by 5 mismatches in the hom,
  // then followed by 20 matches in the right side.
  std::string     querySeq1 = "TCTATCACCCATCGTACCACTCACGGGAGCTCTCC";
  SRAlignmentInfo alignmentInfo1;
  alignmentInfo1.leftSize  = 10;
  alignmentInfo1.homSize   = 5;
  alignmentInfo1.rightSize = 20;
  // Here mismatches present only in the hom.
  calculateAlignScore(querySeq1, targetSeq, 9, alignmentInfo1);
  BOOST_REQUIRE_EQUAL(alignmentInfo1.alignScore, 30);

  // Alignment of following query sequence starts at position 9 (0-based)
  // with 5 matches, 3 mismatches, 2 matches in the left side followed by 5 matches
  // in the hom, then followed by 20 matches in the right side.
  std::string     querySeq2 = "TCTATGTTCCTATTAACCACTCACGGGAGCTCTCC";
  SRAlignmentInfo alignmentInfo2;
  alignmentInfo2.leftSize  = 10;
  alignmentInfo2.homSize   = 5;
  alignmentInfo2.rightSize = 20;
  // Here mismatches present only in the left side.
  calculateAlignScore(querySeq2, targetSeq, 9, alignmentInfo2);
  BOOST_REQUIRE_EQUAL(alignmentInfo2.alignScore, 32);

  // Alignment of following query sequence starts at position 9 (0-based)
  // with 10 matches in the left side followed by 5 matches in the hom,
  // then followed by 13 matches,2 mismatches,1 match, 2 mismatches,2 matches
  // in the right side.
  std::string     querySeq3 = "TCTATCACCCTATTAACCACTCACGGGATGTGACC";
  SRAlignmentInfo alignmentInfo3;
  alignmentInfo3.leftSize  = 10;
  alignmentInfo3.homSize   = 5;
  alignmentInfo3.rightSize = 20;
  // Here mismatches present only in the right side.
  calculateAlignScore(querySeq3, targetSeq, 9, alignmentInfo3);
  BOOST_REQUIRE_EQUAL(alignmentInfo3.alignScore, 31);

  // Alignment of following query sequence starts at position 9 (0-based)
  // with 3 matches, 1 mismatches,1match, 1 mismatches,4 matches in the left side
  // followed by 5 mismatches in the hom, then followed by 14 matches, 1 mismatches,
  // 5 matches in the right side.
  std::string     querySeq4 = "TCTGTTACCCATCGTACCACTCACGGGAGTTCTCC";
  SRAlignmentInfo alignmentInfo4;
  alignmentInfo4.leftSize  = 10;
  alignmentInfo4.homSize   = 5;
  alignmentInfo4.rightSize = 20;
  // Here mismatches present in the three sides.
  calculateAlignScore(querySeq4, targetSeq, 9, alignmentInfo4);
  BOOST_REQUIRE_EQUAL(alignmentInfo4.alignScore, 27);

  SRAlignmentInfo alignmentInfo5;
  // Alignment of following query sequence starts at position 9 (0-based)
  // with 35 matches.
  std::string querySeq5    = "TCTATCACCCTATTAACCACTCACGGGAGCTCTCC";
  alignmentInfo5.leftSize  = 10;
  alignmentInfo5.homSize   = 5;
  alignmentInfo5.rightSize = 20;
  // Here no mismatch in any of the location
  calculateAlignScore(querySeq5, targetSeq, 9, alignmentInfo5);
  BOOST_REQUIRE_EQUAL(alignmentInfo5.alignScore, 35);
}

// Test the probability likelihood of a read alignment to the target sequence
// Let's say snp prior probability is sp, then for each base
// 1) log((sp-e/3) + (1-e)*sp) + log(1/3), if the query base is not matched with the aligned target base
// 2) log((sp-e/3) + (1-e)*sp) + log(1/4), if the query base or the target base is N
// 3) log(1 - ((sp-e/3) + (1-e)*sp)), if the query base is matched with the aligned target base
//  where e is the probability of the target base being an error
// so the likelyhood of a read is the sum of above values.
BOOST_AUTO_TEST_CASE(test_getLnLhood)
{
  std::string targetSeq =
      "GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGT"
      "ATTTTCGTCTGGGGGGTGTGCACGCGATAGCATTGCGAGACGCTGGA";

  // Alignment of following query sequence starts at position 9 (0-based)
  // with 10 matches followed by 5 mismatches then followed by 20 matches.
  std::string querySeq1 = "TCTATCACCCATCGTACCACTCACGGGAGCTCTCC";
  unsigned    querySize(querySeq1.size());
  // set the base qualities
  std::unique_ptr<uint8_t[]> qual(new uint8_t[querySize]);
  for (unsigned i(0); i < querySize; ++i) {
    qual[i] = 30;
  }
  CallOptionsShared optionsShared;
  qscore_snp        qscoreSnp(optionsShared.snpPrior);  // snp prior proability
  // Likelyhood score will be calculated from location 8
  // to location 49
  known_pos_range2 range(8, 50);

  static const float lnOneThird(std::log(1 / 3.f));
  float              lnlhood1(0);
  for (unsigned i(0); i < querySize; i++) {
    // 1st 10 bases and last 20 bases are match
    if (i < 10 || i > 14)
      lnlhood1 += qscoreSnp.qphred_to_ln_comp_error_prob(30);
    else  // base-10 to base-14 total 5 mismatch bases
      lnlhood1 += qscoreSnp.qphred_to_ln_error_prob(30) + lnOneThird;
  }
  static const float eps = 0.00000001f;
  BOOST_REQUIRE_CLOSE(
      getLnLhood(querySeq1, qscoreSnp, qual.get(), targetSeq, 9, range, false, 0.f), lnlhood1, eps);

  // Alignment of following query sequence starts at position 9 (0-based)
  // with 10 matches followed by 5 mismatches followed by 1 mismatch with N then
  // followed by 19 matches.
  std::string querySeq2 = "TCTATCACCCATCGTNCCACTCACGGGAGCTCTCC";
  float       lnlhood2(0);
  float       lnRandomBase(-std::log(4.f));
  for (unsigned i(0); i < querySize; i++) {
    // 1st 10 bases and last 19 bases are match
    if (i < 10 || i > 15)
      lnlhood2 += qscoreSnp.qphred_to_ln_comp_error_prob(30);
    else if (i == 15)  // base-15 is N
      lnlhood2 += lnRandomBase;
    else  // base-10 to base-14, total 5 mismatch bases
      lnlhood2 += qscoreSnp.qphred_to_ln_error_prob(30) + lnOneThird;
  }
  BOOST_REQUIRE_CLOSE(
      getLnLhood(querySeq2, qscoreSnp, qual.get(), targetSeq, 9, range, false, 0.f), lnlhood2, eps);
}

// Test the alignment information
// Given a reference contig and a bam record, this test should check the
// alignment information with respect to reference.
// Also this test verifies the result of likelyhood which is mentioned in
// the test_getLnLhood.
BOOST_AUTO_TEST_CASE(test_getRefAlignment)
{
  reference_contig_segment refSeq;
  refSeq.seq() =
      "GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGT"
      "ATTTTCGTCTGGGGGGTGTGCACGCGATAGCATTGCGAGACGCTGGA";
  bam_record bamRecord1;
  buildTestBamRecord(bamRecord1, 0, 9, 0, 100, 35, 15, "35M", "TCTATCACCCATCGTACCACTCACGGGAGCTCTCC", 200);

  // alignment information will be calculated using the following range
  //                                | <-range1->|
  //                                |  [18,24)  |
  // Alignment->    --------------------------------------------->
  //                | left size(10) |hom size(5)| right size(20) |
  known_pos_range2   range1(18, 24);
  CallOptionsShared  optionsShared;
  qscore_snp         qscoreSnp(optionsShared.snpPrior);  // snp prior probability
  SRAlignmentInfo    srAlignmentInfo1;
  static const float lnOneThird(std::log(1 / 3.f));

  float lnlhood1(0);
  for (unsigned i(0); i < 35; i++) {
    // 1st 10 bases and last 20 bases are match
    if (i < 10 || i > 14)
      lnlhood1 += qscoreSnp.qphred_to_ln_comp_error_prob(30);
    else  // base-10 to base-14, total 5 mismatch bases
      lnlhood1 += qscoreSnp.qphred_to_ln_error_prob(30) + lnOneThird;
  }

  static const float eps = 0.00000001f;
  getRefAlignment(bamRecord1, refSeq, range1, qscoreSnp, srAlignmentInfo1);

  // Based on above schematic diagram and value of range1, following values
  // are expected.
  BOOST_REQUIRE_EQUAL(srAlignmentInfo1.leftSize, 10);
  // There is no mismatches in left side
  BOOST_REQUIRE_EQUAL(srAlignmentInfo1.leftMismatches, 0);
  // size of the right side is 20
  BOOST_REQUIRE_EQUAL(srAlignmentInfo1.rightSize, 20);
  // There is no mismatches in right side
  BOOST_REQUIRE_EQUAL(srAlignmentInfo1.rightMismatches, 0);
  // Hom size is 5
  BOOST_REQUIRE_EQUAL(srAlignmentInfo1.homSize, 5);
  // There are 5 mismatches in hom size
  BOOST_REQUIRE_EQUAL(srAlignmentInfo1.homMismatches, 5);
  BOOST_REQUIRE_CLOSE(srAlignmentInfo1.alignLnLhood, lnlhood1, eps);

  // alignment information will be calculated using the following range
  //                                | <-range1->|
  //                                |  [19,24)  |
  // Alignment->    --------------------------------------------->
  //                | left size(11) |hom size(4)| right size(20) |
  known_pos_range2 range2(19, 24);
  SRAlignmentInfo  srAlignmentInfo2;
  getRefAlignment(bamRecord1, refSeq, range2, qscoreSnp, srAlignmentInfo2);
  BOOST_REQUIRE_EQUAL(srAlignmentInfo2.leftSize, 11);
  BOOST_REQUIRE_EQUAL(srAlignmentInfo2.leftMismatches, 1);
  BOOST_REQUIRE_EQUAL(srAlignmentInfo2.rightSize, 20);
  BOOST_REQUIRE_EQUAL(srAlignmentInfo2.rightMismatches, 0);
  BOOST_REQUIRE_EQUAL(srAlignmentInfo2.homSize, 4);
  BOOST_REQUIRE_EQUAL(srAlignmentInfo2.homMismatches, 4);
  BOOST_REQUIRE_CLOSE(srAlignmentInfo2.alignLnLhood, lnlhood1, eps);

  // alignment information will be calculated using the following range
  //                                | <-range1->|
  //                                |  [19,23)  |
  // Alignment->    --------------------------------------------->
  //                | left size(11) |hom size(3)| right size(21) |
  //
  // alignLnLhood remains unchanged as it is calculated on whole read, but
  // here only homology range is changing. So all the sizes will have
  // effect.
  known_pos_range2 range3(19, 23);
  SRAlignmentInfo  srAlignmentInfo3;
  getRefAlignment(bamRecord1, refSeq, range3, qscoreSnp, srAlignmentInfo3);
  BOOST_REQUIRE_EQUAL(srAlignmentInfo3.leftSize, 11);
  BOOST_REQUIRE_EQUAL(srAlignmentInfo3.leftMismatches, 1);
  BOOST_REQUIRE_EQUAL(srAlignmentInfo3.rightSize, 21);
  BOOST_REQUIRE_EQUAL(srAlignmentInfo3.rightMismatches, 1);
  BOOST_REQUIRE_EQUAL(srAlignmentInfo3.homSize, 3);
  BOOST_REQUIRE_EQUAL(srAlignmentInfo3.homMismatches, 3);
  BOOST_REQUIRE_CLOSE(srAlignmentInfo3.alignLnLhood, lnlhood1, eps);

  // Alignment of following query sequence starts at position 9 (0-based)
  // with 10 matches followed by 5 mismatches followed by 1 mismatch with N then
  // followed by 19 matches followed by 5 insertion.
  bam_record bamRecord2;
  buildTestBamRecord(
      bamRecord2, 0, 9, 0, 100, 40, 15, "35M5I", "TCTATCACCCATCGTNCCACTCACGGGAGCTCTCCAGCTA", 200);
  float lnlhood2(0);
  float lnRandomBase(-std::log(4.f));
  for (unsigned i(0); i < 35; i++) {
    // 1st 10 bases and last 19 bases are match
    if (i < 10 || i > 15)
      lnlhood2 += qscoreSnp.qphred_to_ln_comp_error_prob(30);
    else if (i == 15)  // base-15 is N
      lnlhood2 += lnRandomBase;
    else  // base-10 to base-14, total 5 mismatch bases
      lnlhood2 += qscoreSnp.qphred_to_ln_error_prob(30) + lnOneThird;
  }

  known_pos_range2 range4(18, 24);
  SRAlignmentInfo  srAlignmentInfo4;
  getRefAlignment(bamRecord2, refSeq, range4, qscoreSnp, srAlignmentInfo4);
  BOOST_REQUIRE_EQUAL(srAlignmentInfo4.leftSize, 10);
  // There is no mismatches in left side
  BOOST_REQUIRE_EQUAL(srAlignmentInfo4.leftMismatches, 0);
  // size of the right side is 20
  BOOST_REQUIRE_EQUAL(srAlignmentInfo4.rightSize, 20);
  // There is 1 mismatche(base N in the read) in right side
  BOOST_REQUIRE_EQUAL(srAlignmentInfo4.rightMismatches, 1);
  // Hom size is 5
  BOOST_REQUIRE_EQUAL(srAlignmentInfo4.homSize, 5);
  // There are 5 mismatches in hom size
  BOOST_REQUIRE_EQUAL(srAlignmentInfo4.homMismatches, 5);
  BOOST_REQUIRE_CLOSE(srAlignmentInfo4.alignLnLhood, lnlhood2, eps);
}

// Test the exceptions
// Following points need to be tested:
// 1. Query sequence size should not be equal or more than target sequence size
// 2. Scan start should not be be more than scan end where scan start and scan end
// constructs an overlapping window of break point region
BOOST_AUTO_TEST_CASE(test_Exception)
{
  SRAlignmentInfo alignmentInfo1;
  std::string     targetSeq1 = "GATCACAGGTCTAT";
  std::string     querySeq1  = "TCTATCACCCATCGTACCACTCACGGGAGCTCTCC";
  unsigned        querySize(querySeq1.size());
  // set the base qualities
  std::unique_ptr<uint8_t[]> qual(new uint8_t[querySize]);
  for (unsigned i(0); i < querySize; ++i) {
    qual[i] = 30;
  }
  CallOptionsShared optionsShared;
  qscore_snp        qscoreSnp(optionsShared.snpPrior);  // snp probability
  known_pos_range2  range1(8, 50);
  // target sequence is more than query sequence, it will throw an exception
  BOOST_CHECK_THROW(
      splitReadAligner(2, querySeq1, qscoreSnp, qual.get(), targetSeq1, range1, alignmentInfo1),
      illumina::common::GeneralException);

  SRAlignmentInfo alignmentInfo2;
  std::string     targetSeq2 =
      "GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGT"
      "ATTTTCGTCTGGGGGGTGTGCACGCGATAGCATTGCGAGACGCTGGA";
  known_pos_range2 range2(1000, 1050);
  // scan start = max(0, range1.beg() - querysize + 2 ) = 967 and
  // scan end = max(0, min(range1.end(), targetSize - querySize)) = 67
  // scanEnd < scanStart.
  BOOST_CHECK_THROW(
      splitReadAligner(2, querySeq1, qscoreSnp, qual.get(), targetSeq2, range2, alignmentInfo2),
      illumina::common::GeneralException);

  // all good means case-1 and case-2 are satisfied. All the metrics are calculated
  // as described in test_getRefAlignment
  known_pos_range2 range3(8, 50);
  splitReadAligner(2, querySeq1, qscoreSnp, qual.get(), targetSeq2, range3, alignmentInfo2);
  BOOST_REQUIRE_EQUAL(alignmentInfo2.leftSize, 0);
  BOOST_REQUIRE_EQUAL(alignmentInfo2.leftMismatches, 0);
  BOOST_REQUIRE_EQUAL(alignmentInfo2.rightSize, 34);
  BOOST_REQUIRE_EQUAL(alignmentInfo2.rightMismatches, 27);
  BOOST_REQUIRE_EQUAL(alignmentInfo2.homSize, 1);
  BOOST_REQUIRE_EQUAL(alignmentInfo2.homMismatches, 1);
  BOOST_REQUIRE_EQUAL(alignmentInfo2.alignScore, 7);
}

BOOST_AUTO_TEST_SUITE_END()
