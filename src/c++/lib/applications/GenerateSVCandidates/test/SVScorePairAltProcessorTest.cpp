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
#include "SVScorePairAltProcessor.hh"
#include "SVScorePairAltProcessor.cpp"

/// TestSVScorerAltProcessor is a friend of SVScorer. So that it can access private
/// method of SVScorePairAltProcessor
struct TestSVScorerAltProcessor
{
    bool alignShadowRead(SVScorePairAltProcessor& altProcessor, bam_record& bamRecord, int& altTemplateSize)
    {
        return altProcessor.alignShadowRead(bamRecord, altTemplateSize);
    }
};

BOOST_AUTO_TEST_SUITE( SVScorePairAltProcessor_test_suite )

// Test chromosome label from bam header
BOOST_AUTO_TEST_CASE( test_getChromName )
{
    const bam_header_info bamHeader(buildTestBamHeader());
    BOOST_REQUIRE_EQUAL(getChromName(bamHeader, 0), "chrFoo");
    // -1 means unknown chromosome
    BOOST_REQUIRE_EQUAL(getChromName(bamHeader, -1), "UNKNOWN");
}

// Test whether a record should be skipped or not.
// For a large insert SV (insert sequence length >= 100)
// 1. Ignore unpaired reads
// 2. Ignore a read which is unmapped as well as mate is also unmapped.
// For a small insert SV (insert sequence length < 100)
// 1. Ignore a read which is mapped but mate is unmapped
// 2. Ignore a read which is unmapped
// 3. Ignore a pair which is not innie pair
// 4. Do not ignore a read which is mapped, its mate is also mapped and pair is innie pair
BOOST_AUTO_TEST_CASE( test_isSkipRecord )
{
    const bam_header_info bamHeader(buildTestBamHeader());
    const std::vector<bool> initIsAlignmentTumor;
    std::unique_ptr<SVLocusScanner> scanner(buildTestSVLocusScanner(bamHeader));
    const PairOptions options1(false);
    SVCandidate candidate;
    // large insert SV. Insert sequence size is greater than 100.
    candidate.insertSeq = "GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGT"
                          "ATTTTCGTCTGGGGGGTGTGCACGCGATAGCATTGCGAGACGCTGGA";
    candidate.bp1.interval.range = known_pos_range2(100, 220);
    candidate.bp2.interval.range = known_pos_range2(250,370);
    candidate.bp1.state = SVBreakendState::RIGHT_OPEN;
    candidate.bp2.state = SVBreakendState::LEFT_OPEN;
    candidate.setPrecise();
    candidate.assemblyAlignIndex = 0;
    SVEvidence evidence;
    evidence.samples.resize(1);

    ReadScannerOptions scannerOptions;
    SVRefinerOptions refinerOptions;
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
    jumpAlignmentResultType.align1 = alignment1;
    jumpAlignmentResultType.align2 = alignment2;
    jumpAlignmentResultType.jumpRange = 2;
    candidateAssemblyData.spanningAlignments.push_back(jumpAlignmentResultType);
    candidateAssemblyData.isSpanning = true;
    candidateAssemblyData.extendedContigs.push_back("GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGT"
                                                    "ATTTTCGTCTGGGGGGTGTGCACGCGATAGCATTGCGAGACGCTGGA");
    std::shared_ptr<SVScorePairAltProcessor> processor1(new SVScorePairAltProcessor(bamHeader, scannerOptions,
                                                        refinerOptions, initIsAlignmentTumor,
                                                        scanner.operator*(), options1, candidateAssemblyData,
                                                        candidate, true, evidence));
    std::string querySeq1 = "TCTATCACCCATTTTACCACTCACGGGAGCTCTCC";
    // bamRecord1 is not paired read. This read should be skipped.
    bam_record bamRecord1;
    buildTestBamRecord(bamRecord1);
    bamRecord1.toggle_is_paired();
    bamRecord1.set_qname("bamRecord1");
    BOOST_REQUIRE(processor1.operator*().isSkipRecord(bamRecord1));

    // bamRecord2 and its mate both are unmapped. This read should be skipped.
    bam_record bamRecord2;
    buildTestBamRecord(bamRecord1);
    bamRecord2.toggle_is_unmapped();
    bamRecord2.toggle_is_mate_unmapped();
    bamRecord2.set_qname("bamRecord2");
    BOOST_REQUIRE(processor1.operator*().isSkipRecord(bamRecord2));

    // Small insert SV (insert sequence size is less than 100)
    candidate.insertSeq = "GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTC";
    std::shared_ptr<SVScorePairAltProcessor> processor2(new SVScorePairAltProcessor(bamHeader, scannerOptions,
                                                        refinerOptions, initIsAlignmentTumor,
                                                        scanner.operator*(), options1, candidateAssemblyData,
                                                        candidate, true, evidence));
    // BamRecord3 is mapped and its mate is unmapped
    bam_record bamRecord3;
    buildTestBamRecord(bamRecord3, 0, 9, -1, -1, 0, 15, "35M", querySeq1);
    bamRecord3.toggle_is_first();
    bamRecord3.set_qname("bamRecord3");
    bamRecord3.toggle_is_mate_unmapped();
    BOOST_REQUIRE(processor2.get()->isSkipRecord(bamRecord3));

    // BamRecord4 is unmapped
    bam_record bamRecord4;
    buildTestBamRecord(bamRecord4);
    bamRecord4.toggle_is_unmapped();
    bamRecord4.set_qname("bamRecord4");
    BOOST_REQUIRE(processor2.get()->isSkipRecord(bamRecord4));

    // BamRecord5 and its mate are translocated pair that means it is
    // not an innie pair.
    bam_record bamRecord5;
    buildTestBamRecord(bamRecord5, 0, 9, 1, 25, 0, 15, "35M", querySeq1);
    bamRecord5.toggle_is_first();
    bamRecord5.set_qname("bamRecord5");
    BOOST_REQUIRE(processor2.get()->isSkipRecord(bamRecord5));

    // BamRecord6 and its mate are mapping in a proper innie pair
    bam_record bamRecord6;
    buildTestBamRecord(bamRecord6, 0, 9, 0, 100, 150, 15, "35M", querySeq1);
    bamRecord6.toggle_is_first();
    bamRecord6.set_qname("bamRecord6");
    BOOST_REQUIRE(!processor2.get()->isSkipRecord(bamRecord6));
}

// Test the following cases for large insertion SV (insert sequence length >= 100):
// 1. If a bamrecord is mapped but it's mate is unmapped, it will not add any fragment support unless
//    next record of the mapped record is its unmapped mate and it is properly shadow aligned.
// 2. Direction of mapped read should match with the breakpoint direction.
// 3. Check for shadow alignment
// 4. Fragment length should be between 50 and 125.
// 5. minimum of (BP1 center pos - fragment start + 1) and (fragment end - BP2 center pos)
//    should be greater than minimum fragment threshold which is 50
BOOST_AUTO_TEST_CASE( test_processClearedRecord_LargeInsertion )
{
    const bam_header_info bamHeader(buildTestBamHeader());
    const std::vector<bool> initIsAlignmentTumor;
    std::unique_ptr<SVLocusScanner> scanner(buildTestSVLocusScanner(bamHeader));
    const PairOptions options1(false);
    SVCandidate candidate1;
    // large insertion (size = 102 > 100)
    candidate1.insertSeq = "GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGT"
                          "ATTTTCGTCTGGGGGGTGTGCACGCGATAGCATTGCGAGACGCTGGA";
    candidate1.bp1.interval.range = known_pos_range2(100, 220);
    candidate1.bp2.interval.range = known_pos_range2(250, 300);
    candidate1.bp1.state = SVBreakendState::RIGHT_OPEN;
    candidate1.bp2.state = SVBreakendState::LEFT_OPEN;
    candidate1.setPrecise();
    candidate1.assemblyAlignIndex = 0;
    SVEvidence evidence;
    evidence.samples.resize(1);

    ReadScannerOptions scannerOptions;
    SVRefinerOptions refinerOptions;
    SVCandidateAssemblyData candidateAssemblyData;
    candidateAssemblyData.bestAlignmentIndex = 0;
    Alignment alignment1;
    alignment1.beginPos = 30;
    std::string testCigar1("35=");
    cigar_to_apath(testCigar1.c_str(), alignment1.apath);
    Alignment alignment2;
    alignment2.beginPos = 30;
    std::string testCigar2("35=");
    cigar_to_apath(testCigar2.c_str(), alignment2.apath);
    SVCandidateAssemblyData::JumpAlignmentResultType jumpAlignmentResultType;
    jumpAlignmentResultType.align1 = alignment1;
    jumpAlignmentResultType.align2 = alignment2;
    jumpAlignmentResultType.jumpRange = 2;
    candidateAssemblyData.spanningAlignments.push_back(jumpAlignmentResultType);
    candidateAssemblyData.isSpanning = true;
    candidateAssemblyData.extendedContigs.push_back("GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGG"
                                                    "TGATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGT"
                                                    "TGATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGT"
                                                    "ATTTTCGTCTGGGGGGTGTGCACGCGATAGCATTGCGAGACGCTGGA");
    // Processor for Breakpoint-1
    std::shared_ptr<SVScorePairAltProcessor> processor1(
            new SVScorePairAltProcessor(bamHeader, scannerOptions, refinerOptions, initIsAlignmentTumor,
                                        scanner.operator*(), options1, candidateAssemblyData,
                                        candidate1, true, evidence));
    SVId id;
    SupportFragments suppFrags;

    std::string querySeq1 = "TCTATCACCCATTTTACCACTCACGGGAGCTCTCC";
    // Case-1 is designed here.
    // bamRecord1 is paired, but mate is unmapped. Till this point mate is not processed. As a
    // result, this fragment is not supporting allele on BP1.
    bam_record bamRecord1;
    buildTestBamRecord(bamRecord1);
    bamRecord1.set_qname("Read1");
    bamRecord1.toggle_is_mate_unmapped();
    processor1.get()->nextBamIndex(0);
    processor1.get()->processClearedRecord(id, bamRecord1, suppFrags);
    BOOST_REQUIRE(!evidence.getSampleEvidence(0)[bamRecord1.qname()].ref.bp1.isFragmentSupport);
    // case-2 is designed here. bamRecord2 is a shadow read of bamRecord1 means it's mate is aligned
    // and this is unmapped. Anchored read should be in forward strand for BP1 as BP1 is right_open. As a result,
    // this fragment is not supporting allele on BP1.
    bam_record bamRecord2;
    buildTestBamRecord(bamRecord2);
    bamRecord2.set_qname("Read1");
    bamRecord2.toggle_is_unmapped();
    processor1.get()->processClearedRecord(id, bamRecord2, suppFrags);
    BOOST_REQUIRE(!evidence.getSampleEvidence(0)[bamRecord2.qname()].ref.bp1.isFragmentSupport);

    // bamRecord2 is satisfied for shadow read criteria. But it is not satisfied all the
    // conditions of shadow alignment. Shadow alignment test cases have been described in test_alignShadowRead.
    bamRecord1.set_qname("Read2");
    bamRecord2.set_qname("Read2");
    bamRecord2.toggle_is_mate_fwd_strand();
    processor1.get()->processClearedRecord(id, bamRecord1, suppFrags);
    processor1.get()->processClearedRecord(id, bamRecord2, suppFrags);
    BOOST_REQUIRE(!evidence.getSampleEvidence(0)[bamRecord2.qname()].ref.bp1.isFragmentSupport);

    // Case-4 and Case-5 have been designed here
    SVCandidate candidate2;
    candidate2.insertSeq = "GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGT"
                           "ATTTTCGTCTGGGGGGTGTGCACGCGATAGCATTGCGAGACGCTGGA";
    // BP1 center pos = 54
    candidate2.bp1.interval.range = known_pos_range2(50, 60);
    // BP2 center pos = 70
    candidate2.bp2.interval.range = known_pos_range2(70, 71);
    candidate2.bp1.state = SVBreakendState::RIGHT_OPEN;
    candidate2.bp2.state = SVBreakendState::LEFT_OPEN;
    candidate2.setPrecise();
    candidate2.assemblyAlignIndex = 0;
    // Fragment length = 124
    // minimum of (BP1 center pos - fragment start + 1) and (fragment end - BP2 center pos)
    // = min(54-4, 128-70) = 50 >= min fragment size threshold 50.
    bam_record bamRecord3;
    buildTestBamRecord(bamRecord3, 0, 4, 0, 100, 124, 15, "35M", querySeq1);
    bamRecord3.set_qname("bamRecord3");
    std::shared_ptr<SVScorePairAltProcessor> processor2(
            new SVScorePairAltProcessor(bamHeader, scannerOptions, refinerOptions, initIsAlignmentTumor,
                                        scanner.operator*(), options1, candidateAssemblyData, candidate2,
                                        true, evidence));
    processor2.get()->nextBamIndex(0);
    processor2.get()->processClearedRecord(id, bamRecord3, suppFrags);
    BOOST_REQUIRE(evidence.getSampleEvidence(0)[bamRecord3.qname()].alt.bp1.isFragmentSupport);
    BOOST_REQUIRE(evidence.getSampleEvidence(0)[bamRecord3.qname()].ref.bp1.isFragmentSupport);
}

// Test whether a bam read satisfies the following criteria for small insertion (length less than 100)
// 1. Start of bam read should overlap with the search range where search range is calculated as:
//           search range start = BP center pos - (max fragment size of the sample - min fragment support)
//           search range end = BP center pos + (max fragment size of the sample - min fragment support) + 1
// 2. Alt Fragment length should be greater than or equal to min fragment length of the sample
// 3. Alt Fragment length should be less than or equal to max fragment length of the sample
//        Where Alt Fragment length = fragment length - alt_shift
//                                  alt_shift = BP2 center pos - BP1 center pos - insert sequence size
// 4. minimum of (BP1 center pos - fragment start + 1) and (fragment end - BP2 center pos)
//    should be greater than minimum fragment threshold which is 50
BOOST_AUTO_TEST_CASE( test_processClearedRecord_smallInsertion )
{
    const bam_header_info bamHeader(buildTestBamHeader());
    const std::vector<bool> bamFileInfo = {false};
    std::unique_ptr<SVLocusScanner> scanner(buildTestSVLocusScanner(bamHeader));
    const PairOptions options1(false);
    SVCandidate candidate;
    candidate.insertSeq = "GATCACAGGTCTATCACCCTATTAACCACTC";
    // BP1 center pos = 159
    candidate.bp1.interval.range = known_pos_range2(100, 220);
    // BP2 center pos = 274
    candidate.bp2.interval.range = known_pos_range2(250, 300);
    candidate.bp1.state = SVBreakendState::RIGHT_OPEN;
    candidate.bp2.state = SVBreakendState::LEFT_OPEN;
    candidate.setPrecise();
    candidate.assemblyAlignIndex = 0;
    SVEvidence evidence;
    evidence.samples.resize(1);

    ReadScannerOptions scannerOptions;
    SVRefinerOptions refinerOptions;
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
    jumpAlignmentResultType.align1 = alignment1;
    jumpAlignmentResultType.align2 = alignment2;
    jumpAlignmentResultType.jumpRange = 2;
    candidateAssemblyData.spanningAlignments.push_back(jumpAlignmentResultType);
    candidateAssemblyData.isSpanning = true;
    candidateAssemblyData.extendedContigs.push_back("GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGT"
                                                    "ATTTTCGTCTGGGGGGTGTGCACGCGATAGCATTGCGAGACGCTGGA");
    std::shared_ptr<SVScorePairAltProcessor> processor(
            new SVScorePairAltProcessor(bamHeader, scannerOptions, refinerOptions, bamFileInfo,
                                        scanner.operator*(), options1, candidateAssemblyData, candidate,
                                        true, evidence));
    // As we are interested in BP1,
    // search range start = 159 - (125-50) = 84
    // search range end = 159 + (125-50) + 1 = 235
    processor.get()->nextBamIndex(0);
    SVId id;
    SupportFragments suppFrags;

    std::string querySeq1 = "TCTATCACCCATTTTACCACTCACGGGAGCTCTCC";
    // bam read start = 10. It is not overlapping with search range [84, 235).
    // As a result of this, the fragment is not supporting allele on BP1.
    bam_record bamRecord0;
    buildTestBamRecord(bamRecord0, 0, 10, 0, 125, 200, 15, "35M", querySeq1);
    bamRecord0.set_qname("bamRecord0");
    processor.get()->processClearedRecord(id, bamRecord0, suppFrags);
    // Bam read start = 109. It is overlapping with the search range.
    // Here altFragment length = 100 - (274-159-31) = 16 which is less
    // than minimum fragment length(50) of the sample
    // As a result of this, the fragment is not supporting allele on BP1.
    bam_record bamRecord1;
    buildTestBamRecord(bamRecord1, 0, 109, 0, 175, 100, 15, "35M", querySeq1);
    bamRecord1.set_qname("bamRecord1");
    processor.get()->processClearedRecord(id, bamRecord1, suppFrags);
    BOOST_REQUIRE(!evidence.getSampleEvidence(0)[bamRecord1.qname()].alt.bp1.isFragmentSupport);
    BOOST_REQUIRE(!evidence.getSampleEvidence(0)[bamRecord1.qname()].ref.bp1.isFragmentSupport);

    // Bam read start = 109. It is overlapping with the search range.
    // Here altFragment length = 300 - (274-159-31) = 216 which is greater
    // than maximum fragment length(125) of the
    // sample. As a result of this, the fragment is not supporting allele on BP1.
    bam_record bamRecord2;
    buildTestBamRecord(bamRecord2, 0, 109, 0, 175, 300, 15, "35M", querySeq1);
    bamRecord2.set_qname("bamRecord2");
    processor.get()->processClearedRecord(id, bamRecord2, suppFrags);
    BOOST_REQUIRE(!evidence.getSampleEvidence(0)[bamRecord2.qname()].alt.bp1.isFragmentSupport);
    BOOST_REQUIRE(!evidence.getSampleEvidence(0)[bamRecord2.qname()].ref.bp1.isFragmentSupport);

    // Here case-4 is not satisfied.
    // Here min(159-109+1, 308-274) = 34 which is less than 50
    // As a result of this, the fragment is not supporting allele on BP1.
    bam_record bamRecord3;
    buildTestBamRecord(bamRecord3, 0, 109, 0, 175, 200, 15, "35M", querySeq1);
    bamRecord3.set_qname("bamRecord3");
    processor.get()->processClearedRecord(id, bamRecord3, suppFrags);
    BOOST_REQUIRE(!evidence.getSampleEvidence(0)[bamRecord3.qname()].alt.bp1.isFragmentSupport);
    BOOST_REQUIRE(!evidence.getSampleEvidence(0)[bamRecord3.qname()].ref.bp1.isFragmentSupport);

    candidate.bp1.interval.range = known_pos_range2(80, 81);
    candidate.bp2.interval.range = known_pos_range2(160, 161);
    candidate.insertSeq = "";
    std::shared_ptr<SVScorePairAltProcessor> processor1(new SVScorePairAltProcessor(bamHeader,
                                                        scannerOptions, refinerOptions, bamFileInfo,
                                                        scanner.operator*(), options1, candidateAssemblyData,
                                                        candidate, true, evidence));
    // As we are interested in BP1,
    // search range start = 80 - (125-50) = 5
    // search range end = 80 + (125-50) + 1 = 156
    processor1.get()->nextBamIndex(0);
    id.localId = "INS_1";
    // All the above 4 points are satisfied,
    // bamRecord4 is overlapping with the [5,156).
    // altFragment length = 200 - (160-80-0) = 120 which is greater than 50 and
    // less than 125.
    // Also minimum of (BP1 center pos - fragment start + 1) and (fragment end - BP2 center pos)
    //              = min(80-10+1, 210-160) = 50 >= min fragment length 50.
    // So, the fragment is supporting allele on BP1.
    // It will add spanning pair(PR) support for this segment
    bam_record bamRecord4;
    buildTestBamRecord(bamRecord4, 0, 10, 0, 125, 200, 15, "35M", querySeq1);
    bamRecord4.set_qname("bamRecord4");
    processor1.get()->processClearedRecord(id, bamRecord4, suppFrags);
    BOOST_REQUIRE(evidence.getSampleEvidence(0)[bamRecord4.qname()].alt.bp1.isFragmentSupport);
    BOOST_REQUIRE_EQUAL(suppFrags.getSupportFragment(bamRecord4).read1.SVs.size(), 1);
    BOOST_REQUIRE_EQUAL(*(suppFrags.getSupportFragment(bamRecord4).read1.SVs["INS_1"].begin()), "PR");
    BOOST_REQUIRE_EQUAL(suppFrags.getSupportFragment(bamRecord4).read2.SVs.size(), 1);
    BOOST_REQUIRE_EQUAL(*(suppFrags.getSupportFragment(bamRecord4).read2.SVs["INS_1"].begin()), "PR");
    BOOST_REQUIRE(evidence.getSampleEvidence(0)[bamRecord4.qname()].ref.bp1.isFragmentSupport);
}

// Shadow Read - A unmapped paired-end or mate-pair read whose mate is mapped,
//              the unmapped read should have RNAME and POS identical to its mate.
// Shadow alignment means mate rescue where aligner tries to align unmapped read near by to its mate.
// Test the following cases:
// 1. Minimum clipped read length should be at least 40
// 2. Minimum alignment score is at least 0.85.
// 3. shadow alignment in forward strand
// 4. shadow alignment in reverse strand.
// 5. For empty bam record it should throw an exception
BOOST_AUTO_TEST_CASE( test_alignShadowRead )
{
    const bam_header_info bamHeader(buildTestBamHeader());
    const std::vector<bool> initIsAlignmentTumor;
    std::unique_ptr<SVLocusScanner> scanner(buildTestSVLocusScanner(bamHeader));
    const PairOptions options1(false);
    // This is a hack to reduce min fragment size threshold.
    pos_t  *minS((pos_t *)(&options1.minFragSupport));
    *minS = 30;

    SVCandidate candidate;
    candidate.insertSeq = "GATCACAGGTCTATCACCCTATTAACCACTC";
    candidate.bp1.interval.range = known_pos_range2(40, 41);
    candidate.bp2.interval.range = known_pos_range2(50, 51);
    candidate.bp1.state = SVBreakendState::RIGHT_OPEN;
    candidate.bp2.state = SVBreakendState::LEFT_OPEN;
    candidate.setPrecise();
    candidate.assemblyAlignIndex = 0;
    SVEvidence evidence;
    evidence.samples.resize(1);

    // Creating few objects to construct SVScorePairAltProcessor object
    ReadScannerOptions scannerOptions;
    SVRefinerOptions refinerOptions;
    SVCandidateAssemblyData candidateAssemblyData;
    candidateAssemblyData.bestAlignmentIndex = 0;
    Alignment alignment1;
    alignment1.beginPos = 30;
    std::string testCigar1("50=");
    cigar_to_apath(testCigar1.c_str(), alignment1.apath);
    Alignment alignment2;
    alignment2.beginPos = 30;
    std::string testCigar2("50=");
    cigar_to_apath(testCigar2.c_str(), alignment2.apath);
    SVCandidateAssemblyData::JumpAlignmentResultType jumpAlignmentResultType;
    jumpAlignmentResultType.align1 = alignment1;
    jumpAlignmentResultType.align2 = alignment2;
    jumpAlignmentResultType.jumpRange = 2;
    candidateAssemblyData.spanningAlignments.push_back(jumpAlignmentResultType);
    candidateAssemblyData.isSpanning = true;
    candidateAssemblyData.extendedContigs.push_back("GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGT"
                                                    "ATTTTCGTCTGGGGGGTGTGCACGCGATAGCATTGCGAGACGCTGGA");
    std::shared_ptr<SVScorePairAltProcessor> processor(new SVScorePairAltProcessor(
                                                       bamHeader, scannerOptions, refinerOptions, initIsAlignmentTumor,
                                                       scanner.operator*(), options1, candidateAssemblyData,
                                                       candidate, true, evidence));
    processor.get()->nextBamIndex(0);
    std::string querySeq1 = "TCTATCACCCATTTTACCACTCACGGGAGCTCTCC";

    // bamRecord1 is a shadow read.
    // Designed case-1 where read length is 35 which is less than threshold(40)
    // So, shadow alignment is not possible here.
    bam_record bamRecord1;
    buildTestBamRecord(bamRecord1, 0, 4, 0, 125, 0, 15, "", querySeq1);
    TestSVScorerAltProcessor altProcessor;
    int altTemplateSize;
    BOOST_REQUIRE(!altProcessor.alignShadowRead(processor.operator*(), bamRecord1, altTemplateSize));

    // bamRecord2 is a shadow read.
    // Designed case-2. Here alignment score is 0.13 which less than threshold (0.85)
    // It is not a good shadow read.
    bam_record bamRecord2;
    buildTestBamRecord(bamRecord2, 0, 4, 0, 125, 0, 15, "", "TCTATCACCCATTTTACCACTCACGGGAGCTCTCCCATTTTACCACTCAC");
    BOOST_REQUIRE(!altProcessor.alignShadowRead(processor.operator*(), bamRecord2, altTemplateSize));

    // This is a shadow read.
    // Designed case-2 where alignment score is 1.
    // Perfect shadow alignment. Read sequence is matching with the portion of
    // alt contig sequence.
    bam_record bamRecord3;
    buildTestBamRecord(bamRecord3, 0, 30, 0, 30, 0, 0, "", "GGTATTTTCGTCTGGGGGGTGTGCACGCGATAGCATTGCGAGACGCTGGA");
    BOOST_REQUIRE(altProcessor.alignShadowRead(processor.operator*(), bamRecord3, altTemplateSize));


    candidate.bp1.interval.range = known_pos_range2(200, 201);
    candidate.bp2.interval.range = known_pos_range2(300, 301);
    candidateAssemblyData.extendedContigs.clear();
    candidateAssemblyData.extendedContigs.push_back("GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGG"
                                                    "TGATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGT"
                                                    "TGATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGT"
                                                    "ATTTTCGTCTGGGGGGTGTGCACGCGATAGCATTGCGAGACGCTGGA");
    std::shared_ptr<SVScorePairAltProcessor> processor1(new SVScorePairAltProcessor(
                                                        bamHeader, scannerOptions, refinerOptions, initIsAlignmentTumor,
                                                        scanner.operator*(), options1, candidateAssemblyData,
                                                        candidate, true, evidence));
    // Anchored read in forward strand. So need to reverse complement of the shadow read.
    // Here 50 bases of the read matches with last 50 bases of the contig sequence in the reverse
    // complement fashion. Perfect shadow alignment in reverse strand.
    bam_record bamRecord4;
    buildTestBamRecord(bamRecord4, 0, 4, 0, 125, 0, 0, "", "TCCAGCGTCTCGCAATGCTATCGCGTGCACACCCCCCAGACGAAAATACC");
    bamRecord4.toggle_is_mate_fwd_strand();
    BOOST_REQUIRE(altProcessor.alignShadowRead(processor1.operator*(), bamRecord4, altTemplateSize));

    // For empty record it should throw an exception
    bam_record bamRecord5;
    buildTestBamRecord(bamRecord5);
    std::string emptyString = "";
    std::unique_ptr<uint8_t[]> qual(new uint8_t[0]);
    edit_bam_read_and_quality(emptyString.c_str(), qual.get(), *(bamRecord5.get_data()));
    BOOST_CHECK_THROW(altProcessor.alignShadowRead(processor.operator*(), bamRecord5, altTemplateSize),
                      illumina::common::GeneralException);

    // shadow bam read
    bam_record bamRecord6;
    buildTestBamRecord(bamRecord6, 0, 4, 0, 125, 0, 0, "", "AACAGCGTCTCGCAATGCTATCGCGTGCACACCCCCCAGACGAAAATATT");
    bamRecord6.toggle_is_mate_fwd_strand();
    candidate.insertSeq = "";
    // These insertions haven't been assembled all the way through
    // So there will be some clipping in the alignment. Clipped length
    // should be greater than 40.
    candidate.isUnknownSizeInsertion = true;
    std::shared_ptr<SVScorePairAltProcessor> processor3(new SVScorePairAltProcessor(
                                                        bamHeader, scannerOptions, refinerOptions, initIsAlignmentTumor,
                                                        scanner.operator*(), options1, candidateAssemblyData,
                                                        candidate, true, evidence));
    // Clipped length should be greater than 40
    // Clipping is happening at the end.
    // Alignment is 1=49S. clipped size = 50 - 49 = 1 < 40
    BOOST_REQUIRE(!altProcessor.alignShadowRead(processor3.operator*(), bamRecord6, altTemplateSize));

    candidate.bp1.interval.range = known_pos_range2(40, 41);
    candidate.bp2.interval.range = known_pos_range2(50, 51);
    std::shared_ptr<SVScorePairAltProcessor> processor4(new SVScorePairAltProcessor(
                                                        bamHeader, scannerOptions, refinerOptions, initIsAlignmentTumor,
                                                        scanner.operator*(), options1, candidateAssemblyData,
                                                        candidate, true, evidence));
    // shadow bam read
    bam_record bamRecord7;
    buildTestBamRecord(bamRecord7, 0, 30, 0, 30, 0, 0, "", "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCGGA");
    // Clipped length should be greater than 40
    // Clipping is happening at the beginning.
    // Alignment is 49S1=. clipped size = 50 - 49 = 1 < 40
    BOOST_REQUIRE(!altProcessor.alignShadowRead(processor4.operator*(), bamRecord7, altTemplateSize));
}

BOOST_AUTO_TEST_SUITE_END()
