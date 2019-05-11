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

#include "boost/test/unit_test.hpp"
#include "manta/BamStreamerUtils.hpp"
#include "test/testAlignmentDataUtil.hpp"
#include "test/testFileMakers.hpp"

#include "SVScorerSplit.cpp"

BOOST_AUTO_TEST_SUITE(SVScorerSplit_test_suite)

// Create Temporary bam streams of a bam file which contains
// two bam record.
struct BamStream {
  BamStream()
  {
    const bam_header_info bamHeader(buildTestBamHeader());

    std::string querySeq1 = "TCTATCACCCATTTTACCACTCACGGGAGCTCTCC";
    bam_record  bamRecord1;
    buildTestBamRecord(bamRecord1, 0, 9, 0, 90, 35, 15, "35M", querySeq1, 100);
    bamRecord1.toggle_is_first();
    bamRecord1.set_qname("bamRecord1");
    readsToAdd.push_back(bamRecord1);
    buildTestBamFile(bamHeader, readsToAdd, _bamFilename());

    const std::string                          referenceFilename = getTestReferenceFilename();
    std::vector<std::string>                   bamFilenames      = {_bamFilename()};
    std::vector<std::shared_ptr<bam_streamer>> bamStreams;
    openBamStreams(referenceFilename, bamFilenames, bamStreams);
    bamStream = bamStreams[0];
  }

  std::shared_ptr<bam_streamer> bamStream;
  std::vector<bam_record>       readsToAdd;

private:
  const std::string&     _bamFilename() const { return _bamFilenameMaker.getFilename(); }
  const BamFilenameMaker _bamFilenameMaker;
};

BOOST_FIXTURE_TEST_SUITE(SVSupports_test_suite, BamStream)

// Test the root mean square value of evidence reads
// Following cases need to be tested:
// 1. If there are three reads with mapping qualiity m1, m2, m3 respectively,
//    rms = sqrt(m1*m1 + m2*m2 + m3*m3) / 3
// 2. If split read count is zero, rms value will be zero.
BOOST_AUTO_TEST_CASE(test_finishSampleSRData)
{
  // split read count = 0, so rms value is also 0.
  SVSampleInfo sampleInfo;
  finishSampleSRData(sampleInfo);
  BOOST_REQUIRE_EQUAL(sampleInfo.alt.splitReadMapQ, 0);
  BOOST_REQUIRE_EQUAL(sampleInfo.ref.splitReadMapQ, 0);

  // Given three reads with mapping quality 30 each for alt allele.
  sampleInfo.alt.splitReadMapQ  = 30 * 30 + 30 * 30 + 30 * 30;
  sampleInfo.alt.splitReadCount = 3;
  // Given three reads with mapping quality 40 each for ref allele.
  sampleInfo.ref.splitReadMapQ  = 40 * 40 + 40 * 40 + 40 * 40;
  sampleInfo.ref.splitReadCount = 3;
  finishSampleSRData(sampleInfo);
  BOOST_REQUIRE_EQUAL(sampleInfo.alt.splitReadMapQ, 30);
  BOOST_REQUIRE_EQUAL(sampleInfo.ref.splitReadMapQ, 40);
}

// Test the split read count for ref and alt allele.
// Following cases need to be tested:
// 1. A split read is aligned to the alt contig and the ref contig, and
//    the read support alt more than ref.
// 2. Another split read is aligned to the alt contig and the ref contig, and
//    the read support ref more than alt.
// 3. Sum squared mapping quality value.
BOOST_AUTO_TEST_CASE(test_incrementSplitReadEvidence)
{
  // Assume that we got following information about breakpoint-1 when
  // a read aligns to alt allele contig near breakpoint-1
  SRAlignmentInfo altBp1SR1;
  altBp1SR1.isEvidence   = true;
  altBp1SR1.evidence     = 1.2;
  altBp1SR1.alignLnLhood = -30;
  // Assume that we got following information about breakpoint-2 when
  // a read aligns to alt allele contig near breakpoint-2
  SRAlignmentInfo altBp2SR1;
  altBp2SR1.isEvidence      = true;
  altBp2SR1.isTier2Evidence = true;
  altBp2SR1.evidence        = 1.5;
  altBp2SR1.alignLnLhood    = -20;
  // allele-specific evidence info
  SVSampleAlleleInfo altSVSampleAlleleInfo;
  // track all support data from an individual read in a fragment specific to an
  // individual breakend of a alt allele
  SVFragmentEvidenceAlleleBreakendPerRead altBp1Support1;
  SVFragmentEvidenceAlleleBreakendPerRead altBp2Support1;
  // Assume that we got following information about breakpoint-1 when
  // a read aligns to reference near breakpoint-1
  SRAlignmentInfo refBp1SR1;
  refBp1SR1.isEvidence   = true;
  refBp1SR1.evidence     = 1.4;
  refBp1SR1.alignLnLhood = -40;
  // Assume that we got following information about breakpoint-1 when
  // a read aligns to reference near breakpoint-1
  SRAlignmentInfo refBp2SR1;
  refBp2SR1.isEvidence      = true;
  refBp2SR1.isTier2Evidence = true;
  refBp2SR1.evidence        = 1.4;
  refBp2SR1.alignLnLhood    = -35;
  // ref-specific evidence info
  SVSampleAlleleInfo refSVSampleAlleleInfo;
  // track all support data from an individual read in a fragment specific to an
  // individual breakend of a ref allele
  SVFragmentEvidenceAlleleBreakendPerRead refBp1Support1;
  SVFragmentEvidenceAlleleBreakendPerRead refBp2Support1;
  // Mapping quality = 50
  incrementSplitReadEvidence(
      refBp1SR1,
      refBp2SR1,
      altBp1SR1,
      altBp2SR1,
      50,
      false,
      refSVSampleAlleleInfo,
      altSVSampleAlleleInfo,
      refBp1Support1,
      refBp2Support1,
      altBp1Support1,
      altBp2Support1);
  // snp likelihood score of alt alignment is more than ref alignment,
  // so all the information will be updated for alt alignment.
  BOOST_REQUIRE(altBp1Support1.isSplitSupport);
  BOOST_REQUIRE(!altBp1Support1.isTier2SplitSupport);
  BOOST_REQUIRE_EQUAL(altBp1Support1.splitEvidence, altBp1SR1.evidence);
  BOOST_REQUIRE(altBp2Support1.isSplitSupport);
  BOOST_REQUIRE(altBp2Support1.isTier2SplitSupport);
  BOOST_REQUIRE_EQUAL(altBp2Support1.splitEvidence, altBp2SR1.evidence);
  BOOST_REQUIRE_EQUAL(altSVSampleAlleleInfo.splitReadCount, 1);
  BOOST_REQUIRE_EQUAL(altSVSampleAlleleInfo.splitReadMapQ, 2500);  // sum squared of mapping quality
  BOOST_REQUIRE_EQUAL(
      altSVSampleAlleleInfo.splitReadEvidence, std::max(altBp1SR1.evidence, altBp2SR1.evidence));

  altBp2SR1.alignLnLhood    = -35;
  refBp2SR1.alignLnLhood    = -20;
  refBp1SR1.isTier2Evidence = true;
  refBp2SR1.isTier2Evidence = false;
  // Mapping quality = 20
  incrementSplitReadEvidence(
      refBp1SR1,
      refBp2SR1,
      altBp1SR1,
      altBp2SR1,
      20,
      false,
      refSVSampleAlleleInfo,
      altSVSampleAlleleInfo,
      refBp1Support1,
      refBp2Support1,
      altBp1Support1,
      altBp2Support1);
  // snp likelihood score of ref alignment is more than alt alignment,
  // so all the information will be updated for ref alignment.
  BOOST_REQUIRE(refBp1Support1.isSplitSupport);
  BOOST_REQUIRE(refBp1Support1.isTier2SplitSupport);
  BOOST_REQUIRE_EQUAL(refBp1Support1.splitEvidence, refBp1SR1.evidence);
  BOOST_REQUIRE(refBp2Support1.isSplitSupport);
  BOOST_REQUIRE(!refBp2Support1.isTier2SplitSupport);
  BOOST_REQUIRE_EQUAL(refBp2Support1.splitEvidence, refBp2SR1.evidence);
  BOOST_REQUIRE_EQUAL(refSVSampleAlleleInfo.splitReadCount, 1);
  BOOST_REQUIRE_EQUAL(refSVSampleAlleleInfo.splitReadMapQ, 400);
  BOOST_REQUIRE_EQUAL(
      refSVSampleAlleleInfo.splitReadEvidence, std::max(refBp1SR1.evidence, refBp2SR1.evidence));
}

// The scoreSplitReads api tries to align read on reference regions as well
// alt contigs. Following test verifies the alignment information specially
// the split read evidence count when reference haplotype wins
// over alt haplotype based on alignLnLhood score.
BOOST_AUTO_TEST_CASE(test_AlignmentInfo_When_Reference_Haplotype_Wins)
{
  // Reference sequence at breakpoint-1
  std::string refSeqBp1 =
      "GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGT"
      "ATTTTCGTCTGGGGGGTGTGCACGCGATAGCATTGCGAGACGCTGGA";
  // Alt Contig sequence at breakpoint-1
  std::string altContigSeqBp1 =
      "GATCACAGGTCGATATCCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGT"
      "ATTTTCGTCTGGGGGGTGTGCACGCGATAGCATTGCGAGACGCTGGA";
  // Reference sequence at breakpoint-2
  std::string refSeqBp2 =
      "GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGT"
      "ATTTTCGTCTGGGGGGTGTGCACGCGATAGCATTGCGAGACGCTGGA";
  // Alt Contig sequence at breakpoint-2
  std::string altContigSeqBp2 =
      "GATCACAGGTCGATATCCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGT"
      "ATTTTCGTCTGGGGGGTGTGCACGCGATAGCATTGCGAGACGCTGGA";
  // Options for  snp prior probability
  CallOptionsShared optionsShared;
  // A pair of ids for both ends of a single SV junction
  SVId       id;
  SVBreakend breakend;
  breakend.interval.range = known_pos_range2(24, 25);
  reference_contig_segment reference;
  reference.seq() = refSeqBp1;
  // assembly data for  a specific sv candidate
  SVCandidateAssemblyData svCandidateAssemblyData;
  SVCandidate             candidate;
  // The next four lines are hack to avoid construction of complex svCandidateAssemblyData
  SVAlignmentInfo alignmentInfo(candidate, svCandidateAssemblyData);
  const_cast<std::string&>(alignmentInfo.bp1ReferenceSeq()).assign(refSeqBp1);
  const_cast<std::string&>(alignmentInfo.bp1ContigSeq()).assign(altContigSeqBp1);
  const_cast<std::string&>(alignmentInfo.bp2ReferenceSeq()).assign(refSeqBp2);
  const_cast<std::string&>(alignmentInfo.bp2ContigSeq()).assign(altContigSeqBp2);
  // Contig and reference offset.
  // Gapless alignment starting point will be in between (offset_start - flankScoreSize) and
  // (offset_end + flankScoreSize). In this test case flank score size is 17.
  alignmentInfo.bp1ContigOffset = known_pos_range2(18, 19);
  alignmentInfo.bp1RefOffset    = known_pos_range2(18, 19);
  alignmentInfo.bp2ContigOffset = known_pos_range2(23, 24);
  alignmentInfo.bp2RefOffset    = known_pos_range2(24, 25);

  SVEvidence::evidenceTrack_t evidence;        // this will track all support data for a sv hypothesis
  SVSampleInfo                sample;          // sample specif evidence info
  SVEvidenceWriterSampleData  svSupportFrags;  // set of supporting fragments for an evidence bam
  CallOptionsSharedDeriv      optionsSharedDeriv(optionsShared);

  // Bases on the query sequence and all the breakpoint sequences, following information we will get,
  // AltBP1 Information:
  //         leftSize=10 homSize=1 rightSize=24 leftMismatches=4 homMismatches=1 rightMismatches=1
  //         alignScore=29 alignLnLhood: -48.0602
  // AltBP2 Information:
  //         leftSize=15 homSize=1 rightSize=19 leftMismatches=6 homMismatches=0 rightMismatches=0
  //         alignScore=29 alignLnLhood: -48.0652
  // RefBP1 Information:
  //         leftSize=10 homSize=1 rightSize=24 leftMismatches=1 homMismatches=1 rightMismatches=1
  //         alignScore=32 alignLnLhood: -21.9917
  // RefBP2 Information:
  //         leftSize=16 homSize=1 rightSize=18 leftMismatches=3 homMismatches=0 rightMismatches=0
  //         alignScore=32 alignLnLhood: -22.0037
  // Based on the above information, we can see likelihood score of RefBP1 is more than all, so
  // reference haplotype will be selected and all the information will be updated on ref allele.
  scoreSplitReads(
      optionsSharedDeriv,
      17,
      id,
      breakend,
      alignmentInfo,
      reference,
      true,
      15,
      10,
      0,
      0,
      false,
      evidence,
      bamStream.operator*(),
      sample,
      svSupportFrags);
  // Ref information should be updated
  BOOST_REQUIRE_EQUAL(sample.ref.splitReadCount, 1);
  BOOST_REQUIRE_EQUAL(sample.ref.splitReadMapQ, readsToAdd[0].map_qual() * readsToAdd[0].map_qual());
  BOOST_REQUIRE_EQUAL(sample.alt.splitReadCount, 0);
}

// The scoreSplitReads api tries to align read on reference regions as well
// alt contigs. Following test verifies the alignment information specially
// the split read evidence count when alt haplotype wins
// over reference haplotype based on alignLnLhood score.
BOOST_AUTO_TEST_CASE(test_AlignmentInfo_When_Allele_Haplotype_Wins)
{
  // Reference sequence at breakpoint-1
  std::string refSeqBp1 =
      "GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGT"
      "ATTTTCGTCTGGGGGGTGTGCACGCGATAGCATTGCGAGACGCTGGA";
  // Contig sequence at breakpoint-1
  std::string altContigSeqBp1 =
      "GATCACAGGTCGATATCCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGT"
      "ATTTTCGTCTGGGGGGTGTGCACGCGATAGCATTGCGAGACGCTGGA";
  // Reference sequence at breakpoint-2
  std::string refSeqBp2 =
      "GATCACAGGTCGATATCCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGT"
      "ATTTTCGTCTGGGGGGTGTGCACGCGATAGCATTGCGAGACGCTGGA";
  // Contig sequence at breakpoint-2
  std::string altContigSeqBp2 =
      "GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGT"
      "ATTTTCGTCTGGGGGGTGTGCACGCGATAGCATTGCGAGACGCTGGA";
  // Options for  snp prior probability
  CallOptionsShared optionsShared;
  SVId              id;  // A pair of ids for both ends of a single SV junction
  SVBreakend        breakend;
  breakend.interval.range = known_pos_range2(24, 25);
  reference_contig_segment reference;
  reference.seq() = refSeqBp1;
  SVCandidateAssemblyData svCandidateAssemblyData;
  SVCandidate             candidate;
  // assembly data for  a specific sv candidate
  SVAlignmentInfo alignmentInfo(candidate, svCandidateAssemblyData);
  // The next four lines are hack to avoid construction of complex svCandidateAssemblyData
  const_cast<std::string&>(alignmentInfo.bp1ReferenceSeq()).assign(refSeqBp1);
  const_cast<std::string&>(alignmentInfo.bp1ContigSeq()).assign(altContigSeqBp1);
  const_cast<std::string&>(alignmentInfo.bp2ReferenceSeq()).assign(refSeqBp2);
  const_cast<std::string&>(alignmentInfo.bp2ContigSeq()).assign(altContigSeqBp2);
  // Contig and reference offset
  // Gapless alignment starting point will be in between (offset_start - flankScoreSize) and
  // (offset_end + flankScoreSize). In this test case flank score size is 17.
  alignmentInfo.bp1ContigOffset = known_pos_range2(18, 19);
  alignmentInfo.bp1RefOffset    = known_pos_range2(18, 19);
  alignmentInfo.bp2ContigOffset = known_pos_range2(24, 25);
  alignmentInfo.bp2RefOffset    = known_pos_range2(23, 24);

  SVEvidence::evidenceTrack_t evidence;        // this will track all support data for a sv hypothesis
  SVSampleInfo                sample;          // sample specif evidence info
  SVEvidenceWriterSampleData  svSupportFrags;  // set of supporting fragments for an evidence bam
  CallOptionsSharedDeriv      optionsSharedDeriv(optionsShared);

  // Bases on the query sequence and all the breakpoint sequences, following information we will get,
  // AltBP1 Information:
  //         leftSize=10 homSize=1 rightSize=24 leftMismatches=1 homMismatches=1 rightMismatches=1
  //         alignScore=32 alignLnLhood: -24.0441
  // AltBP2 Information:
  //         leftSize=16 homSize=1 rightSize=18 leftMismatches=3 homMismatches=0 rightMismatches=0
  //         alignScore=32 alignLnLhood: -24.0501
  // RefBP1 Information:
  //         leftSize=10 homSize=1 rightSize=24 leftMismatches=4 homMismatches=1 rightMismatches=1
  //         alignScore=29 alignLnLhood: -43.9273
  // RefBP2 Information:
  //        leftSize=15 homSize=1 rightSize=19 leftMismatches=6 homMismatches=0 rightMismatches=0
  //        alignScore=29  alignLnLhood: -43.9373
  // Based on the above information, we can see likelihood score of AltBP1 is more than all, so
  // alt haplotype will be selected and all the information will be updated on alt allele.
  scoreSplitReads(
      optionsSharedDeriv,
      17,
      id,
      breakend,
      alignmentInfo,
      reference,
      true,
      15,
      10,
      0,
      0,
      false,
      evidence,
      bamStream.operator*(),
      sample,
      svSupportFrags);
  // Alt allele information will be non-zero
  BOOST_REQUIRE_EQUAL(sample.alt.splitReadCount, 1);
  BOOST_REQUIRE_EQUAL(sample.alt.splitReadMapQ, readsToAdd[0].map_qual() * readsToAdd[0].map_qual());
  BOOST_REQUIRE_EQUAL(sample.ref.splitReadCount, 0);
}

// Test the alignment information for RNA read
// Test the following cases:
// 1. RNA alignment info should not be updated based on the likelihood score
// 2. RNA alignement info should be updated for both ref and alt allele whichever
//    satisfies evidence criteria.
BOOST_AUTO_TEST_CASE(test_AlignmentInfo_For_RNA)
{
  // Reference sequence at breakpoint-1
  std::string refSeqBp1 =
      "GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGT"
      "ATTTTCGTCTGGGGGGTGTGCACGCGATAGCATTGCGAGACGCTGGA";
  // Contig sequence at breakpoint-1
  std::string altContigSeqBp1 =
      "GATCACAGGTCGATATCCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGT"
      "ATTTTCGTCTGGGGGGTGTGCACGCGATAGCATTGCGAGACGCTGGA";
  // Reference sequence at breakpoint-2
  std::string refSeqBp2 =
      "GATCACAGGTCTATCACCCATTTTACCACTCACGGGAGCTCTCCATGCATTTGGT"
      "ATTTTCGTCTGGGGGGTGTGCACGCGATAGCATTGCGAGACGCTGGA";
  // Contig sequence at breakpoint-2
  std::string altContigSeqBp2 =
      "GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGT"
      "ATTTTCGTCTGGGGGGTGTGCACGCGATAGCATTGCGAGACGCTGGA";
  // Options for  snp prior probability
  CallOptionsShared optionsShared;
  SVId              id;  // A pair of ids for both ends of a single SV junction
  SVBreakend        breakend;
  breakend.interval.range = known_pos_range2(23, 24);
  reference_contig_segment reference;
  reference.seq() = refSeqBp1;
  // assembly data for  a specific sv candidate
  SVCandidateAssemblyData svCandidateAssemblyData;
  SVCandidate             candidate;
  // The next four lines are hack to avoid construction of complex svCandidateAssemblyData
  SVAlignmentInfo alignmentInfo(candidate, svCandidateAssemblyData);
  const_cast<std::string&>(alignmentInfo.bp1ContigSeq()).assign(altContigSeqBp1);
  const_cast<std::string&>(alignmentInfo.bp2ContigSeq()).assign(altContigSeqBp2);
  // Gapless alignment starting point will be in between (offset_start - flankScoreSize) and
  // (offset_end + flankScoreSize). In this test case flank score size is 17.
  alignmentInfo.bp1ContigOffset = known_pos_range2(18, 19);
  alignmentInfo.bp2ContigOffset = known_pos_range2(24, 25);

  SVEvidence::evidenceTrack_t evidence1;        // this will track all support data for a sv hypothesis
  SVSampleInfo                sample1;          // sample specif evidence info
  SVEvidenceWriterSampleData  svSupportFrags1;  // set of supporting fragments for an evidence bam
  CallOptionsSharedDeriv      optionsSharedDeriv(optionsShared);

  // Bases on the query sequence and all the breakpoint sequences, following information we will get,
  // AltBP1 Information:
  //         leftSize=10 homSize=1 rightSize=24 leftMismatches=1 homMismatches=1 rightMismatches=1
  //         alignScore=32 alignLnLhood: -24.0441
  // AltBP2 Information:
  //         leftSize=16 homSize=1 rightSize=18 leftMismatches=3 homMismatches=0 rightMismatches=0
  //         alignScore=32 alignLnLhood: -24.0501
  // RefBP1 Information:
  //         leftSize=15 homSize=0 rightSize=20 leftMismatches=3 homMismatches=0 rightMismatches=0
  //         alignScore=32 alignLnLhood: -22.0057
  // RefBP2 Information : Not available (As only for Bp1 it will calculate)
  // As in RNA, we will update both alt and ref, so api does not look likeihood score.
  // Based on the above information, we can see only AltBP2 is satisfied all the conditions
  // as mentioned in test_incrementSplitReadEvidence, so all the information
  // will be updated for alt allele.
  scoreSplitReads(
      optionsSharedDeriv,
      17,
      id,
      breakend,
      alignmentInfo,
      reference,
      true,
      15,
      10,
      0,
      0,
      true,
      evidence1,
      bamStream.operator*(),
      sample1,
      svSupportFrags1);
  // Alt allele should be updated
  BOOST_REQUIRE_EQUAL(sample1.alt.splitReadCount, 1);
  BOOST_REQUIRE_EQUAL(sample1.alt.splitReadMapQ, readsToAdd[0].map_qual() * readsToAdd[0].map_qual());
  BOOST_REQUIRE_EQUAL(sample1.ref.splitReadCount, 0);

  SVEvidence::evidenceTrack_t evidence2;
  SVSampleInfo                sample2;
  SVEvidenceWriterSampleData  svSupportFrags2;
  reference.seq()         = refSeqBp2;
  breakend.interval.range = known_pos_range2(24, 25);

  // Bases on the query sequence and all the breakpoint sequences, following information we will get,
  // AltBP1 Information:
  //         leftSize=10 homSize=1 rightSize=24 leftMismatches=1 homMismatches=1 rightMismatches=1
  //         alignScore=32 alignLnLhood: -24.0441
  // AltBP2 Information:
  //         leftSize=16 homSize=1 rightSize=18 leftMismatches=3 homMismatches=0 rightMismatches=0
  //         alignScore=32 alignLnLhood: -24.0501
  // RefBP1 Information: Not Available (As only for Bp2 it will calculate)
  // RefBP2 Information:
  //        leftSize=16 homSize=0 rightSize=19 leftMismatches=0 homMismatches=0 rightMismatches=0
  //        alignScore=35 alignLnLhood: -0.0700233
  // As in RNA, we will update both alt and ref, so api does not look likelihood score.
  // Based on the above information, we can see AltBP2 and RefBP2 are satisfied all the conditions
  // as mentioned in test_incrementSplitReadEvidence, so all the information will
  // be updated for both alt and ref allele.
  scoreSplitReads(
      optionsSharedDeriv,
      17,
      id,
      breakend,
      alignmentInfo,
      reference,
      false,
      15,
      10,
      0,
      0,
      true,
      evidence2,
      bamStream.operator*(),
      sample2,
      svSupportFrags2);
  // Alt and Ref both should be updated
  BOOST_REQUIRE_EQUAL(sample2.ref.splitReadCount, 1);
  BOOST_REQUIRE_EQUAL(sample2.ref.splitReadMapQ, readsToAdd[0].map_qual() * readsToAdd[0].map_qual());
  BOOST_REQUIRE_EQUAL(sample2.alt.splitReadCount, 1);
  BOOST_REQUIRE_EQUAL(sample2.alt.splitReadMapQ, readsToAdd[0].map_qual() * readsToAdd[0].map_qual());
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
