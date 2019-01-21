//
// Manta - Structural Variant and Indel Caller
// Copyright (c) 2013-2018 Illumina, Inc.
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
/// \file
/// \author Atanu Pal
///

#include "boost/test/unit_test.hpp"
#include "test/testAlignmentDataUtil.hh"
#include "test/testSVLocusScanner.hh"

#include "SVScorePairProcessor.hh"
#include "SVScorePairRefProcessor.hh"
#include "SVScorePairRefProcessor.cpp"
#include "SVScorePairProcessor.cpp"

BOOST_AUTO_TEST_SUITE( SVScorePairProcessor_test_suite )

// Test whether a record should be skipped or not.
// Test the following cases:
// 1. Ignore a mapped, paired read whose mate is not mapped.
// 2. Ignore a read which is unmapped.
// 3. Ignore a mapped, paired read whose mate is mapped, but the pair in not innie
// 4. Accept a mapped read of an innie pair (with mapped mate)
BOOST_AUTO_TEST_CASE( test_isSkipRecord )
{
    const std::vector<bool> bamTumorNormalInfo = {false};
    const bam_header_info bamHeader(buildTestBamHeader());
    std::unique_ptr<SVLocusScanner> scanner(buildTestSVLocusScanner(bamHeader));
    const PairOptions options(false);
    SVCandidate candidate;
    candidate.insertSeq = "GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGT"
                          "ATTTTCGTCTGGGGGGTGTGCACGCGATAGCATTGCGAGACGCTGGA";
    SVEvidence evidence;
    SVScorePairRefProcessor processor(bamTumorNormalInfo,
                                      scanner.operator*(), options, candidate, true, evidence);

    std::string querySeq1 = "TCTATCACCCATTTTACCACTCACGGGAGCTCTCC";
    // BamRecord1 is mapped and its mate is unmapped
    bam_record bamRecord1;
    buildTestBamRecord(bamRecord1, 0, 9, 0, 9, 0, 15, "35M", querySeq1);
    bamRecord1.toggle_is_first();
    bamRecord1.set_qname("bamRecord1");
    bamRecord1.toggle_is_mate_unmapped();
    BOOST_REQUIRE(processor.isSkipRecord(bamRecord1));

    // BamRecord2 is unmapped
    bam_record bamRecord2;
    buildTestBamRecord(bamRecord2);
    bamRecord2.toggle_is_unmapped();
    BOOST_REQUIRE(processor.isSkipRecord(bamRecord2));

    // BamRecord3 and its mate are translocated pair that means it is
    // not an innie pair.
    bam_record bamRecord3;
    buildTestBamRecord(bamRecord3, 0, 9, 1, 25, 0, 15, "35M", querySeq1);
    bamRecord3.toggle_is_first();
    bamRecord3.set_qname("bamRecord3");
    BOOST_REQUIRE(processor.isSkipRecord(bamRecord3));

    // BamRecord4 and its mate are mapping in a proper innie pair
    bam_record bamRecord4;
    buildTestBamRecord(bamRecord4, 0, 9, 0, 100, 150, 15, "35M", querySeq1);
    bamRecord4.toggle_is_first();
    bamRecord4.set_qname("bamRecord4");
    BOOST_REQUIRE(!processor.isSkipRecord(bamRecord4));
}

// Test whether a candidate sv contains large insertion or not
// Threshold for large sv is 100.
BOOST_AUTO_TEST_CASE( test_isLargeInsertSV )
{
    const std::vector<bool> bamTumorNormalInfo = {false};
    const bam_header_info bamHeader(buildTestBamHeader());
    std::unique_ptr<SVLocusScanner> scanner(buildTestSVLocusScanner(bamHeader));
    const PairOptions options(false);
    SVCandidate candidate;
    candidate.insertSeq = "GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGT"
                          "ATTTTCGTCTGGGGGGTGTGCACGCGATAGCATTGCGAGACGCTGGA";
    SVEvidence evidence;
    SVScorePairRefProcessor processor(bamTumorNormalInfo,
                                      scanner.operator*(), options, candidate, true, evidence);
    // Size of insert sequence is 102 (>100)
    BOOST_REQUIRE(processor.isLargeInsertSV(candidate));

    candidate.insertSeq = "GATCACAGGTCTATCACCCTATTAACCACTCACG";
    // Size of insert sequence is 34 (<100)
    BOOST_REQUIRE(!processor.isLargeInsertSV(candidate));
}

// Test whether a bam read satisfies core filter or not where core read filter includes
// 1. Ignore PCR duplicate reads
// 2. Ignore filtered reads
// 3. Ignore supplementary reads with no SA tag
// 4. Ignore secondary reads with no SA tag
BOOST_AUTO_TEST_CASE( test_SkipRecordCore)
{
    const std::vector<bool> bamTumorNormalInfo = {false};
    const bam_header_info bamHeader(buildTestBamHeader());
    std::unique_ptr<SVLocusScanner> scanner(buildTestSVLocusScanner(bamHeader));
    const PairOptions options(false);
    SVCandidate candidate;
    candidate.insertSeq = "GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGT"
                          "ATTTTCGTCTGGGGGGTGTGCACGCGATAGCATTGCGAGACGCTGGA";
    SVEvidence evidence;
    SVScorePairRefProcessor processor(bamTumorNormalInfo,
                                      scanner.operator*(), options, candidate, true, evidence);

    // PCR duplicate read
    bam_record bamRecord1;
    buildTestBamRecord(bamRecord1);
    bamRecord1.toggle_is_duplicate();
    BOOST_REQUIRE(processor.isSkipRecordCore(bamRecord1));

    // Filtered Read
    bam_record bamRecord2;
    buildTestBamRecord(bamRecord2);
    bamRecord2.toggle_is_filtered();
    BOOST_REQUIRE(processor.isSkipRecordCore(bamRecord2));

    // Supplementary but SA tag is not there
    bam_record bamRecord3;
    buildTestBamRecord(bamRecord3);
    bamRecord3.toggle_is_supplementary();
    BOOST_REQUIRE(processor.isSkipRecordCore(bamRecord3));

    // Secondary but SA tag is not there
    bam_record bamRecord4;
    buildTestBamRecord(bamRecord4);
    bamRecord4.toggle_is_secondary();
    BOOST_REQUIRE(processor.isSkipRecordCore(bamRecord4));

    // Non strict supplementary read
    // Where non-strict supplementary means
    // 1) Bam flag is saying it is not supplementary alignment
    // 2) Bam flag is saying it is secondary alignment
    // 3) SA tag is present.
    bam_record bamRecord5;
    buildTestBamRecord(bamRecord5);
    bamRecord5.toggle_is_secondary();
    addSupplementaryAlignmentEvidence(bamRecord5);
    BOOST_REQUIRE(processor.isSkipRecordCore(bamRecord5));

    // Proper innie paired read
    std::string querySeq1 = "TCTATCACCCATTTTACCACTCACGGGAGCTCTCC";
    bam_record bamRecord6;
    buildTestBamRecord(bamRecord6, 0, 9, 0, 100, 150, 15, "35M", querySeq1);
    BOOST_REQUIRE(!processor.isSkipRecordCore(bamRecord6));
}

// Test the search range of a sv candidate which is located around centerPos of sv candidate
// So Search range start = Center postion of breakpoint - (max fragment length - min threshold for fragment length),
// Search range end = Center postion of breakpoint + (max fragment length - min threshold for fragment length) + 1
// where min threshold for fragment length = 50 in manta.
BOOST_AUTO_TEST_CASE( test_nextBAMIndex )
{
    const std::vector<bool> bamFileInfo = {false}; // whether normal bam or tumor bam file
    const bam_header_info bamHeader(buildTestBamHeader());
    // Created a locus scanner such that min fragment length = 50 and max fragment length = 125
    std::unique_ptr<SVLocusScanner> scanner(buildTestSVLocusScanner(bamHeader));
    const PairOptions options(false); // false means DNA options
    SVCandidate candidate;
    candidate.insertSeq = "GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGT"
                          "ATTTTCGTCTGGGGGGTGTGCACGCGATAGCATTGCGAGACGCTGGA";
    candidate.bp1.interval.range = known_pos_range2(100, 101);
    candidate.bp2.interval.range = known_pos_range2(150,151);
    SVEvidence evidence;
    // isBP1 = true. So search range is calculated around BP1.
    SVScorePairRefProcessor processor1(bamFileInfo,
                                      scanner.operator*(), options, candidate, true, evidence);
    // center positin of BP1 = 100
    // So range start = 100 - (125-50) = 25
    // range end = 100 + (125-50) + 1 = 176
    GenomeInterval genomeInterval = processor1.nextBamIndex(0);
    BOOST_REQUIRE_EQUAL(genomeInterval.range.begin_pos(), 25);
    BOOST_REQUIRE_EQUAL(genomeInterval.range.end_pos(), 176);

    // isBP1 = false. So search range is calculated around BP2.
    SVScorePairRefProcessor processor2(bamFileInfo,
                                       scanner.operator*(), options, candidate, false,
                                       evidence);
    // center position of BP2 = 150
    // So range start = 150 - (125-50) = 75
    // range end = 150 + (125-50) + 1 = 226
    genomeInterval = processor2.nextBamIndex(0);
    BOOST_REQUIRE_EQUAL(genomeInterval.range.begin_pos(), 75);
    BOOST_REQUIRE_EQUAL(genomeInterval.range.end_pos(), 226);
}

// Test whether a bam read satisfies the following criteria
// 1. Start of bam read should overlap with the search range
// 2. Fragment length should be greater than or equal to min fragment length of the sample
// 3. Fragment length should be less than or equal to max fragment length of the sample
// 4. minimum of (breakpoint center pos - fragment start + 1) and (fragment end - breakpoint center pos)
//    should be greater than minimum fragment threshold which is 50
// 5. For RNA, bam read should be properly paired.
BOOST_AUTO_TEST_CASE( test_processClearedRecord )
{
    const std::vector<bool> bamFileInfo = {false}; // whether normal bam or tumor bam file
    const bam_header_info bamHeader(buildTestBamHeader());
    std::unique_ptr<SVLocusScanner> scanner(buildTestSVLocusScanner(bamHeader));
    const PairOptions options1(false); // false means DNA options
    SVCandidate candidate;
    candidate.insertSeq = "GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGT"
                          "ATTTTCGTCTGGGGGGTGTGCACGCGATAGCATTGCGAGACGCTGGA";
    // BP1 center pos = 159
    candidate.bp1.interval.range = known_pos_range2(100, 220);
    // BP2 center pos = 309
    candidate.bp2.interval.range = known_pos_range2(250, 370);
    SVEvidence evidence;
    evidence.samples.resize(1);
    SVScorePairRefProcessor processor1(bamFileInfo,
                                       scanner.operator*(), options1, candidate, true, evidence);
    // So Search range start = 159 - (125-50) = 84
    // range end = 159 + (125-50) + 1 = 235
    processor1.nextBamIndex(0);
    SVId id;
    SupportFragments suppFrags;

    std::string querySeq1 = "TCTATCACCCATTTTACCACTCACGGGAGCTCTCC";
    // case-1 is designed here.
    // bam read start = 9. It is not overlapping with search range [84, 235).
    // As a result of this, this fragment is not supporting allele on BP1.
    bam_record bamRecord1;
    buildTestBamRecord(bamRecord1, 0, 9, 0, 100, 150, 15, "35M", querySeq1);
    bamRecord1.set_qname("bamRecord1");
    processor1.processClearedRecord(id, bamRecord1, suppFrags);
    BOOST_REQUIRE(!evidence.getSampleEvidence(0)[bamRecord1.qname()].ref.bp1.isFragmentSupport);

    // case-2 is designed here.
    // Bam read start = 109. It is overlapping with the search range. But
    // its fragment length is 49 which is less than minimum fragment length(50) of the sample
    // As a result of this, this fragment is not supporting allele on BP1.
    bam_record bamRecord2;
    buildTestBamRecord(bamRecord2, 0, 109, 0, 125, 49, 15, "35M", querySeq1);
    bamRecord2.set_qname("bamRecord2");
    processor1.processClearedRecord(id, bamRecord2, suppFrags);
    BOOST_REQUIRE(!evidence.getSampleEvidence(0)[bamRecord2.qname()].ref.bp1.isFragmentSupport);

    // Case-3 is designed here.
    // Bam read start = 109. It is overlapping with the search range. But
    // its fragment length is 130 which is greater than maximum fragment length(125) of the
    // sample. As a result of this, this fragment is not supporting allele on BP1.
    bam_record bamRecord3;
    buildTestBamRecord(bamRecord3, 0, 109, 0, 200, 130, 15, "35M", querySeq1);
    bamRecord3.set_qname("bamRecord3");
    processor1.processClearedRecord(id, bamRecord3, suppFrags);
    BOOST_REQUIRE(!evidence.getSampleEvidence(0)[bamRecord3.qname()].ref.bp1.isFragmentSupport);

    // Case-4 is designed here.
    // Here min(159-109+1, 168-159) = 9 which is less than 50. As a result of this, this fragment
    // is not supporting allele on BP1.
    bam_record bamRecord4;
    buildTestBamRecord(bamRecord4, 0, 109, 0, 125, 60, 15, "35M", querySeq1);
    bamRecord4.set_qname("bamRecord4");
    processor1.processClearedRecord(id, bamRecord4, suppFrags);
    BOOST_REQUIRE(!evidence.getSampleEvidence(0)[bamRecord4.qname()].ref.bp1.isFragmentSupport);

    // Here min(159-109+1, 208-159) = 51 which is greater than 50
    // Fragment start = 109 which is overlapping with the search range[84,235).
    // Fragment size = 100 which is greater than 50 and less than 125.
    // All the above 4 points satisfied, this fragment is supporting allele on BP1.
    bam_record bamRecord5;
    buildTestBamRecord(bamRecord5, 0, 109, 0, 200, 100, 15, "35M", querySeq1);
    bamRecord5.set_qname("bamRecord5");
    processor1.processClearedRecord(id, bamRecord5, suppFrags);
    BOOST_REQUIRE(evidence.getSampleEvidence(0)[bamRecord5.qname()].ref.bp1.isFragmentSupport);

    // The following two test cases are for RNA sample
    const PairOptions options2(true);
    SVScorePairRefProcessor processor2(bamFileInfo,
                                       scanner.operator*(), options2,
                                       candidate, true, evidence);
    // For RNA, read should be in proper pair. Following
    // bam read is not in proper pair. As a result of this, this fragment is
    // not supporting allele on BP1.
    bam_record bamRecord6;
    buildTestBamRecord(bamRecord6, 0, 109, 0, 200, 150, 15, "35M", querySeq1);
    bamRecord6.set_qname("bamRecord6");
    processor2.nextBamIndex(0);
    processor2.processClearedRecord(id, bamRecord6, suppFrags);
    BOOST_REQUIRE(!evidence.getSampleEvidence(0)[bamRecord6.qname()].ref.bp1.isFragmentSupport);

    SVScorePairRefProcessor processor3(bamFileInfo,
                                       scanner.operator*(), options2, candidate, true, evidence);
    // All the above points are satisfied for RNA sample.
    // Also fragment is in proper pair.
    // So this fragment is supporting allele on BP1.
    bam_record bamRecord7;
    buildTestBamRecord(bamRecord7, 0, 109, 0, 200, 150, 15, "35M", querySeq1);
    bamRecord7.set_qname("bamRecord7");
    bamRecord7.get_data()->core.flag ^= BAM_FLAG::PROPER_PAIR;
    processor3.nextBamIndex(0);
    processor3.processClearedRecord(id, bamRecord7, suppFrags);
    BOOST_REQUIRE(evidence.getSampleEvidence(0)[bamRecord7.qname()].ref.bp1.isFragmentSupport);

    // count only from the down stream reads. So here mate location
    // is less than this read's start location. Fragment support is
    // calculated based on mate's location and the fragment length.
    // Here min(159-109+1, 208-159) = 51 which is greater than 50
    // Fragment start = 109 which is overlapping with the search range[84,235).
    // Fragment size = 100 which is greater than 50 and less than 125.
    // All the above 4 points satisfied, this fragment is supporting allele on BP1.
    bam_record bamRecord8;
    buildTestBamRecord(bamRecord8, 0, 200, 0, 109, 150, 15, "84M");
    bamRecord8.set_qname("bamRecord8");
    bamRecord8.get_data()->core.flag ^= BAM_FLAG::PROPER_PAIR;
    processor3.processClearedRecord(id, bamRecord8, suppFrags);
    BOOST_REQUIRE(evidence.getSampleEvidence(0)[bamRecord8.qname()].ref.bp1.isFragmentSupport);
}

BOOST_AUTO_TEST_SUITE_END()
