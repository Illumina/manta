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

#include "SVScorePairAltProcessor.cpp"
#include "SVScorePairAltProcessor.hpp"
#include "boost/test/unit_test.hpp"
#include "test/testAlignmentDataUtil.hpp"
#include "test/testSVLocusScanner.hpp"

/// TestSVScorerAltProcessor is a friend of SVScorer. So that it can access private
/// method of SVScorePairAltProcessor
struct TestSVScorerAltProcessor {
  bool alignShadowRead(
      SVScorePairAltProcessor& altProcessor, const bam_record& bamRecord, int& altTemplateSize)
  {
    return altProcessor.alignShadowRead(bamRecord, altTemplateSize);
  }
};

BOOST_AUTO_TEST_SUITE(SVScorePairAltProcessor_test_suite)

// Test chromosome label from bam header
BOOST_AUTO_TEST_CASE(test_getChromName)
{
  const bam_header_info bamHeader(buildTestBamHeader());
  BOOST_REQUIRE_EQUAL(getChromName(bamHeader, 0), "chrFoo");
  // -1 means unknown chromosome
  BOOST_REQUIRE_EQUAL(getChromName(bamHeader, -1), "UNKNOWN");
}

// Test whether a record should be skipped or not.
// For a large insert SV (insert sequence length >= 100)
// 1. Ignore a mapped read which is not paired.
// 2. Ignore a read which is unmapped as well as mate is also unmapped.
// For a small insert SV (insert sequence length < 100)
// 1. Ignore a mapped, paired read whose mate is not mapped.
// 2. Ignore a read which is unmapped.
// 3. Ignore a mapped, paired read whose mate is mapped, but the pair in not innie.
// 4. Accept a mapped read of an innie pair (with mapped mate).
// Only case-1 of small insert SV is implemented here. All other test cases are described in
// test_isSkipRecord in SVScorePairProcessorTest.cpp.
BOOST_AUTO_TEST_CASE(test_isSkipRecord)
{
  const bam_header_info           bamHeader(buildTestBamHeader());
  const std::vector<bool>         tumorNormalInfo = {false};
  std::unique_ptr<SVLocusScanner> scanner(buildTestSVLocusScanner(bamHeader));
  const PairOptions               options1(false);
  SVCandidate                     candidate;
  // The test cases below are for large insertions, which has insert
  // sequence size larger than or equal to 100.
  candidate.insertSeq =
      "GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGT"
      "ATTTTCGTCTGGGGGGTGTGCACGCGATAGCATTGCGAGACGCTGGA";
  candidate.bp1.interval.range = known_pos_range2(100, 220);
  candidate.bp2.interval.range = known_pos_range2(250, 370);
  candidate.bp1.state          = SVBreakendState::RIGHT_OPEN;
  candidate.bp2.state          = SVBreakendState::LEFT_OPEN;
  candidate.setPrecise();
  candidate.assemblyAlignIndex = 0;
  SVEvidence evidence;
  evidence.samples.resize(1);

  ReadScannerOptions      scannerOptions;
  SVRefinerOptions        refinerOptions;
  SVCandidateAssemblyData candidateAssemblyData;
  candidateAssemblyData.bestAlignmentIndex = 0;
  Alignment alignment1;
  alignment1.beginPos = 406;
  std::string testCigar1("94=");
  cigar_to_apath(testCigar1.c_str(), alignment1.apath);
  Alignment alignment2;
  alignment2.beginPos = 510;
  std::string testCigar2("110=75I1D2=");
  cigar_to_apath(testCigar2.c_str(), alignment2.apath);
  SVCandidateAssemblyData::JumpAlignmentResultType jumpAlignmentResultType;
  jumpAlignmentResultType.align1    = alignment1;
  jumpAlignmentResultType.align2    = alignment2;
  jumpAlignmentResultType.jumpRange = 2;
  candidateAssemblyData.spanningAlignments.push_back(jumpAlignmentResultType);
  candidateAssemblyData.isSpanning = true;
  candidateAssemblyData.extendedContigs.push_back(
      "GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGT"
      "ATTTTCGTCTGGGGGGTGTGCACGCGATAGCATTGCGAGACGCTGGA");
  SVScorePairAltProcessor processor1(
      bamHeader,
      scannerOptions,
      refinerOptions,
      tumorNormalInfo,
      scanner.operator*(),
      options1,
      candidateAssemblyData,
      candidate,
      true,
      evidence);
  std::string querySeq1 = "TCTATCACCCATTTTACCACTCACGGGAGCTCTCC";
  // bamRecord1 is not paired read. This read should be skipped.
  bam_record bamRecord1;
  buildTestBamRecord(bamRecord1);
  bamRecord1.toggle_is_paired();
  bamRecord1.set_qname("bamRecord1");
  BOOST_REQUIRE(processor1.isSkipRecord(bamRecord1));

  // bamRecord2 and its mate both are unmapped. This read should be skipped.
  bam_record bamRecord2;
  buildTestBamRecord(bamRecord2);
  bamRecord2.toggle_is_unmapped();
  bamRecord2.toggle_is_mate_unmapped();
  bamRecord2.set_qname("bamRecord2");
  BOOST_REQUIRE(processor1.isSkipRecord(bamRecord2));

  // The test cases below are for small insertions, which has
  // insert sequence size less than 100.
  candidate.insertSeq = "GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTC";
  SVScorePairAltProcessor processor2(
      bamHeader,
      scannerOptions,
      refinerOptions,
      tumorNormalInfo,
      scanner.operator*(),
      options1,
      candidateAssemblyData,
      candidate,
      true,
      evidence);
  // BamRecord3 is mapped and its mate is unmapped
  bam_record bamRecord3;
  buildTestBamRecord(bamRecord3, 0, 9, 0, 9, 35, 15, "35M", querySeq1, 0);
  bamRecord3.toggle_is_first();
  bamRecord3.set_qname("bamRecord3");
  bamRecord3.toggle_is_mate_unmapped();
  BOOST_REQUIRE(processor2.isSkipRecord(bamRecord3));
}

// Test whether a fragment will support a breakpoint for the alt allele if the following conditions are
// satisfied: For large insertion SV (insert sequence length >= 100):
// 1. If a bamrecord is mapped but it's mate is unmapped, it will not add any fragment support unless
//    it is properly shadow aligned (assuming the next record of the mapped record is its unmapped mate)
// 2. Check the position of a shadow read relative to the breakpoint.
//    The mapped mate should be in forward strand (left of insert) for BP1 and
//    reverse strand (right of insert) for BP2.
// 3. Check for shadow alignment
// For all SVs, including large insertions:
// 1. Start of bam read should overlap with the search range where search range is calculated as:
//           search range start = BP center pos - (max fragment size of the sample - min fragment threshold)
//           search range end = BP center pos + (max fragment size of the sample - min fragment threshold) + 1
//           search range of a sv candidate is located around centerPos of sv candidate.
// 2. Alt Fragment length should be greater than or equal to min fragment length of the sample.
// 3. Alt Fragment length should be less than or equal to max fragment length of the sample
//        Where Alt Fragment length = fragment length - alt_shift
//                                  alt_shift = BP2 center pos - BP1 center pos - insert sequence size
// 4. minimum of (BP1 center pos - fragment start + 1) and (fragment end - BP2 center pos)
//    should be greater than minimum fragment threshold which is 50
// If any of the above conditions is not satisfied, then the fragment will not support the breakpoint.
BOOST_AUTO_TEST_CASE(test_processClearedRecord)
{
  const bam_header_info           bamHeader(buildTestBamHeader());
  const std::vector<bool>         tumorNormalInfo = {false};
  std::unique_ptr<SVLocusScanner> scanner(buildTestSVLocusScanner(bamHeader));
  const PairOptions               options1(false);
  SVCandidate                     candidate1;
  // large insertion (size = 102 > 100)
  candidate1.insertSeq =
      "GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGT"
      "ATTTTCGTCTGGGGGGTGTGCACGCGATAGCATTGCGAGACGCTGGA";
  candidate1.bp1.interval.range = known_pos_range2(100, 220);
  candidate1.bp2.interval.range = known_pos_range2(250, 300);
  candidate1.bp1.state          = SVBreakendState::RIGHT_OPEN;
  candidate1.bp2.state          = SVBreakendState::LEFT_OPEN;
  candidate1.setPrecise();
  candidate1.assemblyAlignIndex = 0;
  SVEvidence evidence1;
  evidence1.samples.resize(1);

  ReadScannerOptions      scannerOptions;
  SVRefinerOptions        refinerOptions;
  SVCandidateAssemblyData candidateAssemblyData1;
  candidateAssemblyData1.bestAlignmentIndex = 0;
  Alignment alignment1;
  alignment1.beginPos = 30;
  std::string testCigar1("35=");
  cigar_to_apath(testCigar1.c_str(), alignment1.apath);
  Alignment alignment2;
  alignment2.beginPos = 30;
  std::string testCigar2("35=");
  cigar_to_apath(testCigar2.c_str(), alignment2.apath);
  SVCandidateAssemblyData::JumpAlignmentResultType jumpAlignmentResultType;
  jumpAlignmentResultType.align1    = alignment1;
  jumpAlignmentResultType.align2    = alignment2;
  jumpAlignmentResultType.jumpRange = 2;
  candidateAssemblyData1.spanningAlignments.push_back(jumpAlignmentResultType);
  candidateAssemblyData1.isSpanning = true;
  candidateAssemblyData1.extendedContigs.push_back(
      "GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGG"
      "TGATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGT"
      "TGATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGT"
      "ATTTTCGTCTGGGGGGTGTGCACGCGATAGCATTGCGAGACGCTGGA");
  // Processor for Breakpoint-1
  SVScorePairAltProcessor processor1(
      bamHeader,
      scannerOptions,
      refinerOptions,
      tumorNormalInfo,
      scanner.operator*(),
      options1,
      candidateAssemblyData1,
      candidate1,
      true,
      evidence1);
  SVId                       id;
  SVEvidenceWriterSampleData suppFrags;
  // The test cases below are for large insertions, which has insert sequence size larger than
  // or equal to 100.
  std::string querySeq1 = "TCTATCACCCATTTTACCACTCACGGGAGCTCTCC";
  // Case-1 is designed here.
  // bamRecord1 is mapped, but mate is unmapped. Till this point, the unmapped mate is not processed,
  // assuming it's the next record of the mapped record. As a
  // result, this fragment is not supporting allele on BP1.
  bam_record bamRecord1;
  buildTestBamRecord(bamRecord1);
  bamRecord1.set_qname("Read1");
  bamRecord1.toggle_is_mate_unmapped();
  processor1.nextBamIndex(0);
  processor1.processClearedRecord(id, bamRecord1, suppFrags);
  BOOST_REQUIRE(!evidence1.getSampleEvidence(0)[bamRecord1.qname()].alt.bp1.isFragmentSupport);

  // case-2 is designed here.
  // Check the position of a shadow read relative to the breakpoint.
  // The mapped mate should be in forward strand (left of insert) for BP1 and
  // reverse strand (right of insert) for BP2.
  // Here bamRecord2 is a shadow read (unmapped) while its mate bamRecord1 is mapped.
  // bamrecord1 is in reverse strand, and we are evaluating for BP1.
  // So this fragment is not supporting allele on BP1.
  bam_record bamRecord2;
  buildTestBamRecord(bamRecord2);
  bamRecord2.set_qname("Read1");
  bamRecord2.toggle_is_unmapped();
  processor1.processClearedRecord(id, bamRecord2, suppFrags);
  BOOST_REQUIRE(!evidence1.getSampleEvidence(0)[bamRecord2.qname()].alt.bp1.isFragmentSupport);

  // bamRecord4 meets the basic criteria for shadow read criteria given it is unmapped
  // and its mate bamRecord3 is mapped. However, one of the conditions for shadow alignment
  // is that the clipped read length shall be >= 40. But here it is 35. As a result,
  // this fragment is not supporting allele on BP1.
  // There are many other conditions for shadow alignment which are described in test_alignShadowRead.
  bam_record bamRecord3;
  buildTestBamRecord(bamRecord3);
  bamRecord3.set_qname("Read2");
  bamRecord3.toggle_is_mate_unmapped();

  bam_record bamRecord4;
  buildTestBamRecord(bamRecord4);
  bamRecord4.set_qname("Read2");
  bamRecord4.toggle_is_unmapped();
  bamRecord4.toggle_is_mate_fwd_strand();
  processor1.processClearedRecord(id, bamRecord3, suppFrags);
  processor1.processClearedRecord(id, bamRecord4, suppFrags);
  BOOST_REQUIRE(!evidence1.getSampleEvidence(0)[bamRecord4.qname()].alt.bp1.isFragmentSupport);

  // The test cases below are for all SVs, including large insertions.
  SVCandidate candidate2;
  candidate2.insertSeq = "GATCACAGGTCTATCACCCTATTAACCACTC";
  // BP1 center pos = 159
  candidate2.bp1.interval.range = known_pos_range2(100, 220);
  // BP2 center pos = 274
  candidate2.bp2.interval.range = known_pos_range2(250, 300);
  candidate2.bp1.state          = SVBreakendState::RIGHT_OPEN;
  candidate2.bp2.state          = SVBreakendState::LEFT_OPEN;
  candidate2.setPrecise();
  candidate2.assemblyAlignIndex = 0;
  SVEvidence evidence2;
  evidence2.samples.resize(1);
  SVCandidateAssemblyData candidateAssemblyData2;
  candidateAssemblyData2.bestAlignmentIndex = 0;
  Alignment alignment3;
  alignment3.beginPos = 406;
  std::string testCigar3("94=");
  cigar_to_apath(testCigar3.c_str(), alignment3.apath);
  Alignment alignment4;
  alignment4.beginPos = 510;
  std::string testCigar4("110=75I1D2=");
  cigar_to_apath(testCigar4.c_str(), alignment4.apath);
  jumpAlignmentResultType.align1    = alignment3;
  jumpAlignmentResultType.align2    = alignment4;
  jumpAlignmentResultType.jumpRange = 2;
  candidateAssemblyData2.spanningAlignments.push_back(jumpAlignmentResultType);
  candidateAssemblyData2.isSpanning = true;
  candidateAssemblyData2.extendedContigs.push_back(
      "GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGT"
      "ATTTTCGTCTGGGGGGTGTGCACGCGATAGCATTGCGAGACGCTGGA");
  SVScorePairAltProcessor processor2(
      bamHeader,
      scannerOptions,
      refinerOptions,
      tumorNormalInfo,
      scanner.operator*(),
      options1,
      candidateAssemblyData2,
      candidate2,
      true,
      evidence2);
  // As we are interested in BP1,
  // search range start = 159 - (125-50) = 84
  // search range end = 159 + (125-50) + 1 = 235
  // Test cases for search range are mentioned in test_nextBAMIndex in SVScorePairProcessorTest.cpp
  processor2.nextBamIndex(0);
  SVId                       id2;
  SVEvidenceWriterSampleData suppFrags2;

  // Case-1 is designed here.
  // bam read start = 10. It is not overlapping with search range [84, 235).
  // As a result of this, the fragment is not supporting allele on BP1.
  bam_record bamRecord5;
  buildTestBamRecord(bamRecord5, 0, 10, 0, 125, 35, 15, "35M", querySeq1, 200);
  bamRecord5.set_qname("bamRecord5");
  processor2.processClearedRecord(id2, bamRecord5, suppFrags2);
  BOOST_REQUIRE(!evidence2.getSampleEvidence(0)[bamRecord5.qname()].alt.bp1.isFragmentSupport);

  // Case-2 is designed here.
  // Bam read start = 109. It does overlap with the search range. However,
  // here altFragment length = 100 - (274-159-31) = 16 which is less
  // than minimum fragment length(50) of the sample
  // As a result of this, the fragment is not supporting allele on BP1.
  bam_record bamRecord6;
  buildTestBamRecord(bamRecord6, 0, 109, 0, 175, 35, 15, "35M", querySeq1, 100);
  bamRecord6.set_qname("bamRecord6");
  processor2.processClearedRecord(id, bamRecord6, suppFrags2);
  BOOST_REQUIRE(!evidence2.getSampleEvidence(0)[bamRecord6.qname()].alt.bp1.isFragmentSupport);
  BOOST_REQUIRE(!evidence2.getSampleEvidence(0)[bamRecord6.qname()].ref.bp1.isFragmentSupport);

  // Case-3 is designed here.
  // Bam read start = 109. It does overlap with the search range. However,
  // here altFragment length = 300 - (274-159-31) = 216 which is greater
  // than maximum fragment length(125) of the sample.
  // As a result of this, the fragment is not supporting allele on BP1.
  bam_record bamRecord7;
  buildTestBamRecord(bamRecord7, 0, 109, 0, 175, 35, 15, "35M", querySeq1, 300);
  bamRecord7.set_qname("bamRecord7");
  processor2.processClearedRecord(id2, bamRecord7, suppFrags2);
  BOOST_REQUIRE(!evidence2.getSampleEvidence(0)[bamRecord7.qname()].alt.bp1.isFragmentSupport);
  BOOST_REQUIRE(!evidence2.getSampleEvidence(0)[bamRecord7.qname()].ref.bp1.isFragmentSupport);

  // Here case-4 is not satisfied.
  // Here min(159-109+1, 308-274) = 34 which is less than 50
  // As a result of this, the fragment is not supporting allele on BP1.
  bam_record bamRecord8;
  buildTestBamRecord(bamRecord8, 0, 109, 0, 175, 35, 15, "35M", querySeq1, 200);
  bamRecord8.set_qname("bamRecord8");
  processor2.processClearedRecord(id2, bamRecord8, suppFrags2);
  BOOST_REQUIRE(!evidence2.getSampleEvidence(0)[bamRecord8.qname()].alt.bp1.isFragmentSupport);
  BOOST_REQUIRE(!evidence2.getSampleEvidence(0)[bamRecord8.qname()].ref.bp1.isFragmentSupport);

  candidate2.bp1.interval.range = known_pos_range2(80, 81);
  candidate2.bp2.interval.range = known_pos_range2(160, 161);
  candidate2.insertSeq          = "";
  SVScorePairAltProcessor processor3(
      bamHeader,
      scannerOptions,
      refinerOptions,
      tumorNormalInfo,
      scanner.operator*(),
      options1,
      candidateAssemblyData2,
      candidate2,
      true,
      evidence2);
  // As we are interested in BP1,
  // search range start = 80 - (125-50) = 5
  // search range end = 80 + (125-50) + 1 = 156
  processor3.nextBamIndex(0);
  id2.localId = "INS_1";
  // All the above 4 cases are satisfied,
  // bamRecord9 is overlapping with the [5,156).
  // altFragment length = 200 - (160-80-0) = 120 which is greater than 50 and
  // less than 125.
  // Also minimum of (BP1 center pos - fragment start + 1) and (fragment end - BP2 center pos)
  //              = min(80-10+1, 210-160) = 50 >= min fragment length 50.
  // So, the fragment is supporting allele on BP1.
  // It will add spanning pair(PR) support for this segment
  bam_record bamRecord9;
  buildTestBamRecord(bamRecord9, 0, 10, 0, 125, 35, 15, "35M", querySeq1, 200);
  bamRecord9.set_qname("bamRecord9");
  processor3.processClearedRecord(id2, bamRecord9, suppFrags2);
  BOOST_REQUIRE(evidence2.getSampleEvidence(0)[bamRecord9.qname()].alt.bp1.isFragmentSupport);
  BOOST_REQUIRE_EQUAL(suppFrags2.getSupportFragment(bamRecord9).read1.SVs.size(), 1);
  BOOST_REQUIRE_EQUAL(*(suppFrags2.getSupportFragment(bamRecord9).read1.SVs["INS_1"].begin()), "PR");
  BOOST_REQUIRE_EQUAL(suppFrags2.getSupportFragment(bamRecord9).read2.SVs.size(), 1);
  BOOST_REQUIRE_EQUAL(*(suppFrags2.getSupportFragment(bamRecord9).read2.SVs["INS_1"].begin()), "PR");
  BOOST_REQUIRE(evidence2.getSampleEvidence(0)[bamRecord9.qname()].ref.bp1.isFragmentSupport);
}

// Shadow Read - A unmapped read with its paired mate being mapped,
//              assuming the unmapped read should have RNAME and POS identical to its mate.
// Shadow alignment means mate rescue where the aligner tries to align the unmapped read
// near its anchored (mapped) mate.
// Test the following cases:
// 1. Minimum clipped read length should be at least 40
// 2. Alignment score relative to optimal score should be at least 0.85.
// 3. shadow alignment in reverse strand.
// 4. Bam record with empty sequence should throw an exception
BOOST_AUTO_TEST_CASE(test_alignShadowRead)
{
  const bam_header_info           bamHeader(buildTestBamHeader());
  const std::vector<bool>         tumorNormalInfo = {false};
  std::unique_ptr<SVLocusScanner> scanner(buildTestSVLocusScanner(bamHeader));
  const PairOptions               options1(false);
  SVCandidate                     candidate;
  candidate.insertSeq          = "GATCACAGGTCTATCACCCTATTAACCACTC";
  candidate.bp1.interval.range = known_pos_range2(40, 41);
  candidate.bp2.interval.range = known_pos_range2(50, 51);
  candidate.bp1.state          = SVBreakendState::RIGHT_OPEN;
  candidate.bp2.state          = SVBreakendState::LEFT_OPEN;
  candidate.setPrecise();
  candidate.assemblyAlignIndex = 0;
  SVEvidence evidence;
  evidence.samples.resize(1);

  // Creating few objects to construct SVScorePairAltProcessor object
  ReadScannerOptions      scannerOptions;
  SVRefinerOptions        refinerOptions;
  SVCandidateAssemblyData candidateAssemblyData;
  candidateAssemblyData.bestAlignmentIndex = 0;
  Alignment alignment1;
  alignment1.beginPos = 30;
  std::string testCigar1("50=");
  cigar_to_apath(testCigar1.c_str(), alignment1.apath);
  Alignment alignment2;
  alignment2.beginPos = 100;
  std::string testCigar2("5=10I40=");
  cigar_to_apath(testCigar2.c_str(), alignment2.apath);
  SVCandidateAssemblyData::JumpAlignmentResultType jumpAlignmentResultType;
  jumpAlignmentResultType.align1    = alignment1;
  jumpAlignmentResultType.align2    = alignment2;
  jumpAlignmentResultType.jumpRange = 2;
  candidateAssemblyData.spanningAlignments.push_back(jumpAlignmentResultType);
  candidateAssemblyData.isSpanning = true;
  candidateAssemblyData.extendedContigs.push_back(
      "GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGT"
      "ATTTTCGTCTGGGGGGTGTGCACGCGATAGCATTGCGAGACGCTGGA");
  candidate.bp1.interval.range = known_pos_range2(40, 45);
  candidate.bp2.interval.range = known_pos_range2(54, 55);
  SVScorePairAltProcessor processor1(
      bamHeader,
      scannerOptions,
      refinerOptions,
      tumorNormalInfo,
      scanner.operator*(),
      options1,
      candidateAssemblyData,
      candidate,
      true,
      evidence);
  processor1.nextBamIndex(0);
  std::string querySeq1 = "TCTATCACCCATTTTACCACTCACGGGAGCTCTCC";

  // bamRecord1 is a shadow read.
  // Designed case-1 where read length is 35 which is less than threshold(40)
  // So, shadow alignment is not possible here.
  bam_record bamRecord1;
  buildTestBamRecord(bamRecord1, 0, 125, 0, 125, 35, 15, "", querySeq1);
  TestSVScorerAltProcessor altProcessor;
  int                      altTemplateSize;
  BOOST_REQUIRE(!altProcessor.alignShadowRead(processor1, bamRecord1, altTemplateSize));

  // Designed case-1 when clipping happens at the end of the alignment.
  bam_record bamRecord2;
  buildTestBamRecord(
      bamRecord2, 0, 125, 0, 125, 50, 0, "", "AACAGCGTCTCGCAATGCTATCGCGTGCACACCCCCCAGACGAAAATATT");
  bamRecord2.toggle_is_mate_fwd_strand();
  candidate.insertSeq              = "";
  candidate.isUnknownSizeInsertion = true;
  candidate.bp1.interval.range     = known_pos_range2(200, 201);
  candidate.bp2.interval.range     = known_pos_range2(300, 301);
  candidateAssemblyData.extendedContigs.clear();
  candidateAssemblyData.extendedContigs.push_back(
      "GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGG"
      "TGATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGT"
      "TGATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGT"
      "ATTTTCGTCTGGGGGGTGTGCACGCGATAGCATTGCGAGACGCTGGA");
  SVScorePairAltProcessor processor2(
      bamHeader,
      scannerOptions,
      refinerOptions,
      tumorNormalInfo,
      scanner.operator*(),
      options1,
      candidateAssemblyData,
      candidate,
      true,
      evidence);
  // Clipped length should be greater than 40
  // Clipping is happening at the end.
  // Contig begin offset = 0
  // Contig end offset = alignment1.beginPos + cigar_length -1 = 79
  //
  // clang-format off
    // Contig - ....ACTCACGGGAGCTCATTTGGTGATCACAGGTCTATCACCCTATTA
    // Read -                                                   AACAGCGTCTCGCAATGCTATCGCGTGCACACCCCCCAGACGAAAATATT
  // clang-format on
  //
  // CIGAR - 1=49S
  // So, clipped size = 50 - 49 = 1 < 40
  BOOST_REQUIRE(!altProcessor.alignShadowRead(processor2, bamRecord2, altTemplateSize));

  // Designed case-1 when clipping happens at the beginning of the alignment.
  candidate.bp1.interval.range = known_pos_range2(40, 41);
  candidate.bp2.interval.range = known_pos_range2(50, 51);
  SVScorePairAltProcessor processor3(
      bamHeader,
      scannerOptions,
      refinerOptions,
      tumorNormalInfo,
      scanner.operator*(),
      options1,
      candidateAssemblyData,
      candidate,
      true,
      evidence);
  // shadow bam read
  bam_record bamRecord3;
  buildTestBamRecord(
      bamRecord3, 0, 30, 0, 30, 50, 0, "", "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCGGA");
  // Clipped length should be greater than 40
  // Clipping is happening at the beginning.
  // Contig begin offset = alignment1.beginPos + cigar_length -1 = 79
  // Contig end offset = contig length = 213
  // Contig -                                                  ACCACTCACGGGAGCTCTCCATGC.......
  // Read -   CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCGGA
  // CIGAR - 49S1=
  // So, clipped size = 50 - 49 = 1 < 40
  BOOST_REQUIRE(!altProcessor.alignShadowRead(processor3, bamRecord3, altTemplateSize));

  // bamRecord4 is a shadow read.
  // Designed case-2. Here alignment score relative to optimal score is 0.13 which less than threshold (0.85)
  // It is not a good shadow alignment.
  // Below is the detail of exact match.
  //
  // clang-format off
    // Contig - GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGTATTTTCGTCTGGGGGGTGTGCACGCGATAGCATTGCGAGACGCTGGA
    // Read -            TCTATCACCCATTTTACCACTCACGGGAGCTCTCCCATTTTACCACTCAC
  // clang-format on
  //
  // CIGAR is 10=2X2=1X19=15I1=
  // Match score = 2, mismtach penalty = -8, gap opening penalty = -12, gap extension penalty = -1,
  // clipping penalty = -1.
  // Alignment score is 20 - 16 + 4 - 8 + 38 - 27 + 2 = 13
  // Optimal score is = 50*2 = 100
  // Alignment score relative to optimal score is 13/100 = 0.13 which is less than 0.85
  bam_record bamRecord4;
  buildTestBamRecord(
      bamRecord4, 0, 125, 0, 125, 35, 15, "", "TCTATCACCCATTTTACCACTCACGGGAGCTCTCCCATTTTACCACTCAC");
  BOOST_REQUIRE(!altProcessor.alignShadowRead(processor1, bamRecord4, altTemplateSize));

  // This is a shadow read.
  // Designed case-2 where alignment score should be more than 0.85.
  // Here alignment score is 1.
  // Perfect shadow alignment. Read sequence is matching with the portion of
  // alt contig sequence.
  //
  // clang-format off
    // Below is the detail of exact match.
    // Contig - GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGTATTTTCGTCTGGGGGGTGTGCACGCGATAGCATTGCGAGACGCTGGA
    // Read -                              ACTCACGGGAGCTCTCCATGCATTTGGTATTTTCGTCTGGGGGGTGTGCACGCGATAGCATTGCGAGACGCTGGA
  // clang-format on
  //
  // CIGAR - 75=
  // Alignment score is 75*2 = 150
  // Optimal score is = 75*2 = 150
  // Alignment score relative to optimal score is 150/150 = 1 which is more than 0.85
  candidate.bp1.interval.range     = known_pos_range2(40, 45);
  candidate.bp2.interval.range     = known_pos_range2(54, 55);
  candidate.insertSeq              = "GATCACAGGTCTATCACCCTATTAACCACTC";
  candidate.isUnknownSizeInsertion = false;
  candidateAssemblyData.extendedContigs.clear();
  candidateAssemblyData.extendedContigs.push_back(
      "GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGT"
      "ATTTTCGTCTGGGGGGTGTGCACGCGATAGCATTGCGAGACGCTGGA");
  SVScorePairAltProcessor processor4(
      bamHeader,
      scannerOptions,
      refinerOptions,
      tumorNormalInfo,
      scanner.operator*(),
      options1,
      candidateAssemblyData,
      candidate,
      true,
      evidence);
  bam_record bamRecord5;
  buildTestBamRecord(
      bamRecord5,
      0,
      30,
      0,
      30,
      75,
      0,
      "",
      "ACTCACGGGAGCTCTCCATGCATTTGGTATTTTCGTCTGGGGGGTG"
      "TGCACGCGATAGCATTGCGAGACGCTGGA");
  BOOST_REQUIRE(altProcessor.alignShadowRead(processor4, bamRecord5, altTemplateSize));

  // Designed case-3 where shadow alignment happens in reverse strand.
  // Anchored read in forward strand. So need to reverse complement of the shadow read.
  // Here 50 bases of the read matches with last 50 bases of the contig sequence in the reverse
  // complement fashion. Perfect shadow alignment in reverse strand.
  // Below is the detail of exact match.
  //
  // clang-format off
    // Contig - ....CATTTGGTATTTTCGTCTGGGGGGTGTGCACGCGATAGCATTGCGAGACGCTGGAGATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTC
    // Read -            GAGAGCTCCCGTGAGTGGTTAATAGGGTGATAGACCTGTGATCTCCAGCGTCTCGCAATGCTATCGCGTGCACACCCCCCAGACGAAAATACC
    // Rev complement-   GAGAGCTCCCGTGAGTGGTTAATAGGGTGATAGACCTGTGATCTCCAGCGTCTCGCAATGCTATCGCGTGCACACCCCCCAGACGAAAATACC
  // clang-format on
  //
  // CIGAR - 50=
  // Alignment score is 50*2 = 100
  // Optimal score is = 50*2 = 100
  // Alignment score relative to optimal score is 100/100 = 1 which is more than 0.85
  candidate.bp1.interval.range = known_pos_range2(200, 201);
  candidate.bp2.interval.range = known_pos_range2(300, 301);
  candidateAssemblyData.extendedContigs.clear();
  candidateAssemblyData.extendedContigs.push_back(
      "GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGTGATC"
      "ACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGTTGATCACAG"
      "GTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGTATTTTCGTCTGGGG"
      "GGTGTGCACGCGATAGCATTGCGAGACGCTGGAGATCACAGGTCTATCACCCTATTAACCAC"
      "TCACGGGAGCTCTC");
  SVScorePairAltProcessor processor5(
      bamHeader,
      scannerOptions,
      refinerOptions,
      tumorNormalInfo,
      scanner.operator*(),
      options1,
      candidateAssemblyData,
      candidate,
      true,
      evidence);
  bam_record bamRecord6;
  buildTestBamRecord(
      bamRecord6,
      0,
      125,
      0,
      125,
      50,
      0,
      "",
      "GAGAGCTCCCGTGAGTGGTTAATAGGGTGATAGACCTGTGATCTCCAGCGTCTCG"
      "CAATGCTATCGCGTGCACACCCCCCAGACGAAAATACC");
  bamRecord6.toggle_is_mate_fwd_strand();
  BOOST_REQUIRE(altProcessor.alignShadowRead(processor5, bamRecord6, altTemplateSize));

  // Designed case-4 where empty sequence in bam record is not allowed.
  bam_record bamRecord7;
  buildTestBamRecord(bamRecord7);
  std::string                emptyString = "";
  std::unique_ptr<uint8_t[]> qual(new uint8_t[0]);
  edit_bam_read_and_quality(emptyString.c_str(), qual.get(), *(bamRecord7.get_data()));
  BOOST_CHECK_THROW(
      altProcessor.alignShadowRead(processor1, bamRecord7, altTemplateSize),
      illumina::common::GeneralException);
}

BOOST_AUTO_TEST_SUITE_END()
