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
// 1. Read is mapped but mate is unmapped
// 2. Read is unmapped
// 3. Read pair is not innie pair
// 4. Read is mapped, its mate is also mapped and pair is innie pair
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
    std::shared_ptr<SVScorePairRefProcessor> processor(new SVScorePairRefProcessor(bamTumorNormalInfo,
                                                       scanner.operator*(), options, candidate, true, evidence));

    std::string querySeq1 = "TCTATCACCCATTTTACCACTCACGGGAGCTCTCC";
    // BamRecord1 is mapped and its mate is unmapped
    bam_record bamRecord1;
    buildTestBamRecord(bamRecord1, 0, 9, -1, -1, -1, 15, "35M", querySeq1);
    bamRecord1.toggle_is_first();
    bamRecord1.set_qname("bamRecord1");
    bamRecord1.toggle_is_mate_unmapped();
    BOOST_REQUIRE(processor.get()->isSkipRecord(bamRecord1));

    // BamRecord2 is unmapped
    bam_record bamRecord2;
    buildTestBamRecord(bamRecord2);
    bamRecord2.toggle_is_unmapped();
    BOOST_REQUIRE(processor.get()->isSkipRecord(bamRecord2));

    // BamRecord3 and its mate are translocated pair that means it is
    // an innie pair.
    bam_record bamRecord3;
    buildTestBamRecord(bamRecord3, 0, 9, 1, 25, -1, 15, "35M", querySeq1);
    bamRecord3.toggle_is_first();
    bamRecord3.set_qname("bamRecord3");
    BOOST_REQUIRE(processor.get()->isSkipRecord(bamRecord3));

    // BamRecord4 and its mate are mapping in a proper innie pair
    bam_record bamRecord4;
    buildTestBamRecord(bamRecord4, 0, 9, 0, 100, 150, 15, "35M", querySeq1);
    bamRecord4.toggle_is_first();
    bamRecord4.set_qname("bamRecord4");
    BOOST_REQUIRE(!processor.get()->isSkipRecord(bamRecord4));
}

// Test whether a candidate sv contains large insertion or not
// Threshold for large sv is 100.
BOOST_AUTO_TEST_CASE( test_isLargeInsertSV )
{
    const std::vector<bool> bamTumorNormalInfo;
    const bam_header_info bamHeader(buildTestBamHeader());
    std::unique_ptr<SVLocusScanner> scanner(buildTestSVLocusScanner(bamHeader));
    const PairOptions options(false);
    SVCandidate candidate;
    candidate.insertSeq = "GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGT"
                          "ATTTTCGTCTGGGGGGTGTGCACGCGATAGCATTGCGAGACGCTGGA";
    SVEvidence evidence;
    std::shared_ptr<SVScorePairRefProcessor> processor(new SVScorePairRefProcessor(bamTumorNormalInfo,
                                                       scanner.operator*(), options, candidate, true, evidence));
    // Size of insert sequence is 102 (>100)
    BOOST_REQUIRE(processor.get()->isLargeInsertSV(candidate));
    candidate.insertSeq = "GATCACAGGTCTATCACCCTATTAACCACTCACG";
    // Size of insert sequence is 34 (<100)
    BOOST_REQUIRE(!processor.get()->isLargeInsertSV(candidate));
}

// Test whether a bam read satisfies core filter or not where core read filter includes
// 1. Ignore PCR duplicate reads
// 2. Ignore filtered reads
// 3. Ignore supplementary reads with no SA tag
// 4. Ignore secondary reads with no SA tag
BOOST_AUTO_TEST_CASE( test_SkipRecordCore)
{
    const std::vector<bool> bamTumorNormalInfo;
    const bam_header_info bamHeader(buildTestBamHeader());
    std::unique_ptr<SVLocusScanner> scanner(buildTestSVLocusScanner(bamHeader));
    const PairOptions options(false);
    SVCandidate candidate;
    candidate.insertSeq = "GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGT"
                          "ATTTTCGTCTGGGGGGTGTGCACGCGATAGCATTGCGAGACGCTGGA";
    SVEvidence evidence;
    std::shared_ptr<SVScorePairRefProcessor> processor(new SVScorePairRefProcessor(bamTumorNormalInfo,
                                                       scanner.operator*(), options, candidate, true, evidence));

    // PCR duplicate read
    bam_record bamRecord1;
    buildTestBamRecord(bamRecord1);
    bamRecord1.toggle_is_duplicate();
    BOOST_REQUIRE(processor.get()->isSkipRecordCore(bamRecord1));

    // Filtered Read
    bam_record bamRecord2;
    buildTestBamRecord(bamRecord2);
    bamRecord2.toggle_is_filtered();
    BOOST_REQUIRE(processor.get()->isSkipRecordCore(bamRecord2));

    // Supplementary but SA tag is not there
    bam_record bamRecord3;
    buildTestBamRecord(bamRecord3);
    bamRecord3.toggle_is_supplementary();
    BOOST_REQUIRE(processor.get()->isSkipRecordCore(bamRecord3));

    // Secondary but SA tag is not there
    bam_record bamRecord4;
    buildTestBamRecord(bamRecord4);
    bamRecord4.toggle_is_secondary();
    BOOST_REQUIRE(processor.get()->isSkipRecordCore(bamRecord4));

    // Non strict supplementary read
    bam_record bamRecord5;
    buildTestBamRecord(bamRecord5);
    bamRecord5.toggle_is_secondary();
    addSupplementaryAlignmentEvidence(bamRecord5);
    BOOST_REQUIRE(processor.get()->isSkipRecordCore(bamRecord5));

    // Proper innie paired read
    std::string querySeq1 = "TCTATCACCCATTTTACCACTCACGGGAGCTCTCC";
    bam_record bamRecord6;
    buildTestBamRecord(bamRecord6, 0, 9, 0, 100, 150, 15, "35M", querySeq1);
    BOOST_REQUIRE(!processor.get()->isSkipRecordCore(bamRecord6));
}

// Test the search range of a sv candidate which is to be set around centerPos of sv candidate
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
    // isBP1 = true
    std::shared_ptr<SVScorePairRefProcessor> processor(new SVScorePairRefProcessor(bamFileInfo,
                                                       scanner.operator*(), options, candidate, true, evidence));

    // center positin of BP1 = 100
    // So range start = 100 - (125-50) = 25
    // range end = 100 + (125-50) + 1 = 176
    GenomeInterval genomeInterval = processor.get()->nextBamIndex(0);
    BOOST_REQUIRE_EQUAL(genomeInterval.range.begin_pos(), 25);
    BOOST_REQUIRE_EQUAL(genomeInterval.range.end_pos(), 176);

    // isBP1 = false
    std::shared_ptr<SVScorePairRefProcessor> processor1(new SVScorePairRefProcessor(bamFileInfo,
                                                        scanner.operator*(), options, candidate, false, evidence));
    // center position of BP2 = 150
    // So range start = 150 - (125-50) = 75
    // range end = 150 + (125-50) + 1 = 226
    genomeInterval = processor1.get()->nextBamIndex(0);
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
    std::shared_ptr<SVScorePairRefProcessor> processor(new SVScorePairRefProcessor(bamFileInfo,
                                                       scanner.operator*(), options1, candidate, true, evidence));
    // Search genome interval is [84, 235)
    processor.get()->nextBamIndex(0);
    SVId id;
    SupportFragments suppFrags;

    std::string querySeq1 = "TCTATCACCCATTTTACCACTCACGGGAGCTCTCC";
    // bam read start = 9. It is not overlapping with search range [84, 235).
    // As a result of this, fragment support of this bam read will not be set.
    bam_record bamRecord1;
    buildTestBamRecord(bamRecord1, 0, 9, 0, 100, 150, 15, "35M", querySeq1);
    bamRecord1.set_qname("bamRecord1");
    processor.get()->processClearedRecord(id, bamRecord1, suppFrags);
    BOOST_REQUIRE(!evidence.getSampleEvidence(0)[bamRecord1.qname()].ref.bp1.isFragmentSupport);

    // Bam read start = 109. It is overlapping with the search range. But
    // its fragment length is 49 which is less than minimum fragment length(50) of the sample
    // As a result of this, fragment support of this bam read will not be set.
    bam_record bamRecord2;
    buildTestBamRecord(bamRecord2, 0, 109, 0, 125, 49, 15, "35M", querySeq1);
    bamRecord2.set_qname("bamRecord2");
    processor.get()->processClearedRecord(id, bamRecord2, suppFrags);
    BOOST_REQUIRE(!evidence.getSampleEvidence(0)[bamRecord2.qname()].ref.bp1.isFragmentSupport);

    // Bam read start = 109. It is overlapping with the search range. But
    // its fragment length is 130 which is greater than maximum fragment length(125) of the
    // sample. As a result of this, fragment support of this bam read will not be set.
    bam_record bamRecord3;
    buildTestBamRecord(bamRecord3, 0, 109, 0, 200, 130, 15, "35M", querySeq1);
    bamRecord3.set_qname("bamRecord3");
    processor.get()->processClearedRecord(id, bamRecord3, suppFrags);
    BOOST_REQUIRE(!evidence.getSampleEvidence(0)[bamRecord3.qname()].ref.bp1.isFragmentSupport);

    // Here point-4 is not satisfied. As a result of this, fragment support of this
    // bam read will not be set.
    bam_record bamRecord4;
    buildTestBamRecord(bamRecord4, 0, 109, 0, 125, 60, 15, "35M", querySeq1);
    bamRecord4.set_qname("bamRecord4");
    processor.get()->processClearedRecord(id, bamRecord4, suppFrags);
    BOOST_REQUIRE(!evidence.getSampleEvidence(0)[bamRecord4.qname()].ref.bp1.isFragmentSupport);

    // All the above 4 points satisfied, so we go the fragment support for this read.
    bam_record bamRecord5;
    buildTestBamRecord(bamRecord5, 0, 109, 0, 200, 100, 15, "35M", querySeq1);
    bamRecord5.set_qname("bamRecord5");
    processor.get()->processClearedRecord(id, bamRecord5, suppFrags);
    BOOST_REQUIRE(evidence.getSampleEvidence(0)[bamRecord5.qname()].ref.bp1.isFragmentSupport);

    // The following two test cases are for RNA sample
    const PairOptions options2(true);
    std::shared_ptr<SVScorePairRefProcessor> processor1(new SVScorePairRefProcessor(bamFileInfo,
                                                                                    scanner.operator*(), options2, candidate, true, evidence));
    // For RNA, read should be in proper pair. Following
    // bam read is not in proper pair. As a result of this, fragment support
    // of this bam read will not be set.
    bam_record bamRecord6;
    buildTestBamRecord(bamRecord6, 0, 109, 0, 200, 150, 15, "35M", querySeq1);
    bamRecord6.set_qname("bamRecord6");
    processor1.get()->nextBamIndex(0);
    processor1.get()->processClearedRecord(id, bamRecord6, suppFrags);
    BOOST_REQUIRE(!evidence.getSampleEvidence(0)[bamRecord6.qname()].ref.bp1.isFragmentSupport);

    std::shared_ptr<SVScorePairRefProcessor> processor2(new SVScorePairRefProcessor(bamFileInfo,
                                                        scanner.operator*(), options2, candidate, true, evidence));
    // All the above points are satisfied for RNA sample.
    // So we got evidence support for this bam.
    bam_record bamRecord7;
    buildTestBamRecord(bamRecord7, 0, 109, 0, 200, 150, 15, "35M", querySeq1);
    bamRecord7.set_qname("bamRecord7");
    bamRecord7.get_data()->core.flag ^= BAM_FLAG::PROPER_PAIR;
    processor2.get()->nextBamIndex(0);
    processor2.get()->processClearedRecord(id, bamRecord7, suppFrags);
    BOOST_REQUIRE(evidence.getSampleEvidence(0)[bamRecord7.qname()].ref.bp1.isFragmentSupport);

    // count only from the down stream reads. So here mate location
    // is less than this read's start location. Fragment support is
    // calculated based on mate's location and the fragment length.
    bam_record bamRecord8;
    buildTestBamRecord(bamRecord8, 0, 200, 0, 109, 150, 15, "84M");
    bamRecord8.set_qname("bamRecord8");
    bamRecord8.get_data()->core.flag ^= BAM_FLAG::PROPER_PAIR;
    processor2.get()->processClearedRecord(id, bamRecord8, suppFrags);
    BOOST_REQUIRE(evidence.getSampleEvidence(0)[bamRecord8.qname()].ref.bp1.isFragmentSupport);
}

BOOST_AUTO_TEST_SUITE_END()
