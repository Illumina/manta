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

/// TestSVScorerAltProcessor is a friend of SVScorer. So that can access private
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

// Test whether a bam read satisfies core filter or not where core read filter includes
// 1. test_SkipRecordCore in SVScorePairProcessorTest.cpp
// 2. Ignore unpaired reads
// As case-1 is already tested in SVScorePairProcessorTest.cpp, so only valid read and
// case-2 have been tested here.
BOOST_AUTO_TEST_CASE( test_isSkipRecord )
{
    const bam_header_info bamHeader(buildTestBamHeader());
    const std::vector<bool> initIsAlignmentTumor;
    std::unique_ptr<SVLocusScanner> scanner(buildTestSVLocusScanner(bamHeader));
    const PairOptions options1(false);
    SVCandidate candidate;
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
    std::shared_ptr<SVScorePairAltProcessor> processor(new SVScorePairAltProcessor(bamHeader, scannerOptions, refinerOptions, initIsAlignmentTumor,
                                                                                   scanner.operator*(), options1, candidateAssemblyData, candidate, true, evidence));
    std::string querySeq1 = "TCTATCACCCATTTTACCACTCACGGGAGCTCTCC";
    // bamRecord1 satisfies both case-1 and case-2.
    bam_record bamRecord1;
    buildTestBamRecord(bamRecord1, 0, 9, 0, 125, 150, 15, "35M", querySeq1);
    bamRecord1.toggle_is_first();
    bamRecord1.set_qname("bamRecord1");
    BOOST_REQUIRE(!processor.operator*().isSkipRecord(bamRecord1));

    // bamRecord2 is not paired reads. This read should be skipped.
    bam_record bamRecord2;
    buildTestBamRecord(bamRecord2);
    bamRecord2.toggle_is_paired();
    bamRecord2.set_qname("bamRecord2");
    BOOST_REQUIRE(processor.operator*().isSkipRecord(bamRecord2));
}

// Test the following cases for large insertion SV:
// 1. If a bamrecord is aligned but it's mate is unmapped, it will not add any fragment support unless
//    next record is it's unmapped mapped and it is properly shadow aligned.
// 2. does the shadow occur to the left or right of the insertion? It should occur left for BP1 and It should occur right
//    for BP2.
// 3. Check for shodow alignment
// 4. Fragment length should be between 50 and 125.
// 5. minimum of (breakpoint center pos - fragment start + 1) and (fragment end - breakpoint center pos)
//    should be greater than minimum fragment threshold which is 50
BOOST_AUTO_TEST_CASE( test_LargeInsertion )
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
    // Case-1 is designed hhere.
    // bamRecord1 is paired, but mate is unmapped. Till this point mate is not found.As a
    // result of this, fragment support of this bam read will not be set.
    bam_record bamRecord1;
    buildTestBamRecord(bamRecord1);
    bamRecord1.set_qname("Read1");
    bamRecord1.toggle_is_mate_unmapped();
    processor1.get()->nextBamIndex(0);
    processor1.get()->processClearedRecord(id, bamRecord1, suppFrags);
    BOOST_REQUIRE(!evidence.getSampleEvidence(0)[bamRecord1.qname()].ref.bp1.isFragmentSupport);
    // case-2 is designed here. bamRecord2 is a shadow read of bamRecord1 means it's mate is aligned
    // and this is unmapped. Anchored read should be in forward strand for BP1 computation. As a result of
    // this, fragment support of this bam read will not be set.
    bam_record bamRecord2;
    buildTestBamRecord(bamRecord2);
    bamRecord2.set_qname("Read1");
    bamRecord2.toggle_is_unmapped();
    processor1.get()->processClearedRecord(id, bamRecord2, suppFrags);
    BOOST_REQUIRE(!evidence.getSampleEvidence(0)[bamRecord2.qname()].ref.bp1.isFragmentSupport);

    // bamRecord2 is satisfied for shadow read criteria. But it is satisfied all the
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
    candidate2.bp1.interval.range = known_pos_range2(50, 60);
    candidate2.bp2.interval.range = known_pos_range2(70, 71);
    candidate2.bp1.state = SVBreakendState::RIGHT_OPEN;
    candidate2.bp2.state = SVBreakendState::LEFT_OPEN;
    candidate2.setPrecise();
    candidate2.assemblyAlignIndex = 0;
    bam_record bamRecord3;
    buildTestBamRecord(bamRecord3, 0, 4, 0, 125, 124, 15, "35M", querySeq1);
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

// Test whether a bam read satisfies the following criteria for small insertion
// 1. Start of bam read should overlap with the search range
// 2. Fragment length should be greater than or equal to min fragment length of the sample
// 3. Fragment length should be less than or equal to max fragment length of the sample
// 4. minimum of (breakpoint center pos - fragment start + 1) and (fragment end - breakpoint center pos)
//    should be greater than minimum fragment threshold which is 50
BOOST_AUTO_TEST_CASE( test_smallInsertion )
{
    const bam_header_info bamHeader(buildTestBamHeader());
    const std::vector<bool> bamFileInfo = {false};
    std::unique_ptr<SVLocusScanner> scanner(buildTestSVLocusScanner(bamHeader));
    const PairOptions options1(false);
    SVCandidate candidate;
    candidate.insertSeq = "GATCACAGGTCTATCACCCTATTAACCACTC";
    candidate.bp1.interval.range = known_pos_range2(100, 220);
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
    // GenomeInterval: 0:[84,235)
    processor.get()->nextBamIndex(0);
    SVId id;
    SupportFragments suppFrags;

    std::string querySeq1 = "TCTATCACCCATTTTACCACTCACGGGAGCTCTCC";

    // Bam read start = 109. It is overlapping with the search range. But
    // its fragment length is 49 which is less than minimum fragment length(50) of the sample
    // As a result of this, fragment support of this bam read will not be set.
    bam_record bamRecord1;
    buildTestBamRecord(bamRecord1, 0, 109, 0, 175, 100, 15, "35M", querySeq1);
    bamRecord1.set_qname("bamRecord1");
    processor.get()->processClearedRecord(id, bamRecord1, suppFrags);
    BOOST_REQUIRE(!evidence.getSampleEvidence(0)[bamRecord1.qname()].alt.bp1.isFragmentSupport);
    BOOST_REQUIRE(!evidence.getSampleEvidence(0)[bamRecord1.qname()].ref.bp1.isFragmentSupport);

    // Bam read start = 109. It is overlapping with the search range. But
    // its fragment length is 300 which is greater than maximum fragment length(125) of the
    // sample. As a result of this, fragment support of this bam read will not be set.
    bam_record bamRecord2;
    buildTestBamRecord(bamRecord2, 0, 109, 0, 175, 300, 15, "35M", querySeq1);
    bamRecord2.set_qname("bamRecord2");
    processor.get()->processClearedRecord(id, bamRecord2, suppFrags);
    BOOST_REQUIRE(!evidence.getSampleEvidence(0)[bamRecord2.qname()].alt.bp1.isFragmentSupport);
    BOOST_REQUIRE(!evidence.getSampleEvidence(0)[bamRecord2.qname()].ref.bp1.isFragmentSupport);

    // Here point-4 is not satisfied. As a result of this, fragment support of this
    // bam read will not be set.
    bam_record bamRecord3;
    buildTestBamRecord(bamRecord3, 0, 109, 0, 175, 200, 15, "35M", querySeq1);
    bamRecord3.set_qname("bamRecord3");
    processor.get()->processClearedRecord(id, bamRecord3, suppFrags);
    BOOST_REQUIRE(!evidence.getSampleEvidence(0)[bamRecord3.qname()].alt.bp1.isFragmentSupport);
    BOOST_REQUIRE(!evidence.getSampleEvidence(0)[bamRecord3.qname()].ref.bp1.isFragmentSupport);


    candidate.bp1.interval.range = known_pos_range2(50, 60);
    candidate.bp2.interval.range = known_pos_range2(70, 71);
    std::shared_ptr<SVScorePairAltProcessor> processor1(
            new SVScorePairAltProcessor(bamHeader, scannerOptions, refinerOptions, bamFileInfo,
                                        scanner.operator*(), options1, candidateAssemblyData, candidate,
                                        true, evidence));
    processor1.get()->nextBamIndex(0);
    id.localId = "INS_1";
    // All the above 4 points satisfied, so we will get the fragment support for this read.
    // It will add spanning pair(PR) support for this segment
    bam_record bamRecord4;
    buildTestBamRecord(bamRecord4, 0, 4, 0, 125, 120, 15, "35M", querySeq1);
    bamRecord4.set_qname("bamRecord4");
    processor1.get()->processClearedRecord(id, bamRecord4, suppFrags);
    BOOST_REQUIRE(evidence.getSampleEvidence(0)[bamRecord4.qname()].alt.bp1.isFragmentSupport);
    BOOST_REQUIRE_EQUAL(suppFrags.getSupportFragment(bamRecord4).read1.SVs.size(), 1);
    BOOST_REQUIRE_EQUAL(*(suppFrags.getSupportFragment(bamRecord4).read1.SVs["INS_1"].begin()), "PR");
    BOOST_REQUIRE_EQUAL(suppFrags.getSupportFragment(bamRecord4).read2.SVs.size(), 1);
    BOOST_REQUIRE_EQUAL(*(suppFrags.getSupportFragment(bamRecord4).read2.SVs["INS_1"].begin()), "PR");
    BOOST_REQUIRE(evidence.getSampleEvidence(0)[bamRecord4.qname()].ref.bp1.isFragmentSupport);
}

// Shadow Read-For a unmapped paired-end or mate-pair read whose mate is mapped,
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
    // Designed case-1 wherher read length is 35 which is less than threshold(40)
    // So, shadow alignment is not possible here.
    bam_record bamRecord1;
    buildTestBamRecord(bamRecord1, 0, 4, 0, 125, 0, 15, "", querySeq1);
    TestSVScorerAltProcessor altProcessor;
    int altTemplateSize = 0;
    BOOST_REQUIRE(!altProcessor.alignShadowRead(processor.operator*(), bamRecord1, altTemplateSize));

    // bamRecord2 is a shadow read.
    // Designed case-2. Here alignment score is 0.13 which less than threshold (0.85)
    // It is not a good shadow read.
    bam_record bamRecord2;
    buildTestBamRecord(bamRecord2, 0, 4, 0, 125, 0, 15, "", "TCTATCACCCATTTTACCACTCACGGGAGCTCTCCCATTTTACCACTCAC");
    altTemplateSize = 0;
    BOOST_REQUIRE(!altProcessor.alignShadowRead(processor.operator*(), bamRecord2, altTemplateSize));

    // This is a shadow read.
    // Designed case-2 where alignment score is 1.
    // Perfect shadow alignment. Read sequence is matching with the portion of
    // alt contig sequence.
    bam_record bamRecord3;
    buildTestBamRecord(bamRecord3, 0, 30, 0, 30, 0, 0, "", "GGTATTTTCGTCTGGGGGGTGTGCACGCGATAGCATTGCGAGACGCTGGA");
    altTemplateSize = 0;
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
    // But here 50 bases of the read matches with last 50 bases of the contig sequence in the reverse
    // complement fashion. Perfect shadow alignment in reverse strand.
    bam_record bamRecord4;
    buildTestBamRecord(bamRecord4, 0, 4, 0, 125, 0, 0, "", "TCCAGCGTCTCGCAATGCTATCGCGTGCACACCCCCCAGACGAAAATACC");
    bamRecord4.toggle_is_mate_fwd_strand();
    altTemplateSize = 0;
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
    // these insertions haven't been assembled all the way through
    candidate.isUnknownSizeInsertion = true;
    std::shared_ptr<SVScorePairAltProcessor> processor3(new SVScorePairAltProcessor(
                                                        bamHeader, scannerOptions, refinerOptions, initIsAlignmentTumor,
                                                        scanner.operator*(), options1, candidateAssemblyData,
                                                        candidate, true, evidence));
    altTemplateSize = 0;
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
    altTemplateSize = 0;
    // Alignment is 49S1=. clipped size = 50 - 49 = 1 < 40
    BOOST_REQUIRE(!altProcessor.alignShadowRead(processor4.operator*(), bamRecord7, altTemplateSize));
}

BOOST_AUTO_TEST_SUITE_END()
