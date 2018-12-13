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
//

/// \file
/// \author Atanu Pal
///

#include "boost/test/unit_test.hpp"
#include "test/testAlignmentDataUtil.hh"
#include "test/testSVLocusScanner.hh"
#include "test/testFileMakers.hh"
#include "manta/BamStreamerUtils.hh"

#include "SVScorer.hh"
#include "SVScorerPair.cpp"

/// TestSVScorer is a friend of SVScorer. So that can access private
/// methods of SVScorer.
struct TestSVScorer
{
    void processExistingAltPairInfo(SVScorer &scorer, PairOptions &pairOptions, SVCandidateSetData &candidateSetData,
                                    SVCandidate &candidate, SVId &id, SVEvidence &evidence,
                                    SupportSamples &suppSamples) {
        scorer.processExistingAltPairInfo(pairOptions, candidateSetData, candidate, id, evidence, suppSamples);
    }

    void getSVPairSupport(SVScorer& scorer, SVCandidateSetData& candidateSetData, SVCandidateAssemblyData assemblyData,
                          SVCandidate& candidate, SVId& id, SVEvidence& evidence, SupportSamples& suppSamples)
    {
        scorer.getSVPairSupport(candidateSetData, assemblyData, candidate, id, evidence, suppSamples);
    }
};

BOOST_AUTO_TEST_SUITE( SVScorerPair_test_suite )

// Create Temporary bam streams of a bam file which contains
// three bam records.
struct BamStream
{
    BamStream()
    {
        const bam_header_info bamHeader(buildTestBamHeader());

        std::string querySeq = "TCTATCACCCATTTTACCACTCACGGGAGCTCTCC";
        bam_record bamRecord1;
        buildTestBamRecord(bamRecord1, 0, 9, 0, 70, 125, 15, "35M", querySeq);
        bamRecord1.set_qname("bamRecord1");
        bam_record bamRecord2;
        buildTestBamRecord(bamRecord2, 0, 109, 0, 170, 130, 15, "35M", querySeq);
        bamRecord2.set_qname("bamRecord2");
        bam_record bamRecord3;
        buildTestBamRecord(bamRecord3, 0, 109, 0, 170, 130, 15, "35M", querySeq);
        bamRecord3.set_qname("bamRecord3");
        bam_record bamRecord4;
        buildTestBamRecord(bamRecord4, 0, 109, 0, 175, 100, 15, "35M", querySeq);
        bamRecord4.set_qname("bamRecord4");
        readsToAdd.push_back(bamRecord1);
        readsToAdd.push_back(bamRecord2);
        readsToAdd.push_back(bamRecord3);
        readsToAdd.push_back(bamRecord4);
        bamFileName = _bamFilename();
        buildTestBamFile(bamHeader, readsToAdd, bamFileName);

        const std::string referenceFilename = getTestReferenceFilename();
        std::vector<std::string> bamFilenames = { bamFileName };
        std::vector<std::shared_ptr<bam_streamer>> bamStreams;
        openBamStreams(referenceFilename, bamFilenames, bamStreams);
        bamStream = bamStreams[0];
    }

    std::shared_ptr<bam_streamer> bamStream;
    std::vector<bam_record> readsToAdd;
    std::string bamFileName;
    
private:

    const std::string&
    _bamFilename() const
    {
        return _bamFilenameMaker.getFilename();
    }
    const BamFilenameMaker _bamFilenameMaker;
};

BOOST_FIXTURE_TEST_SUITE( SVScorerPair_test_suite, BamStream )

// getFragInfo api fills read information as accurate as possible based
// on either only mate-1 or mate-2 or both. Following cases need to be tested
// 1. When only mate-1 read available.
// 2. When mate-1 and mate-2 both are available.
BOOST_AUTO_TEST_CASE( test_getFragInfo )
{
    SVCandidateSetSequenceFragment fragment1;
    SpanReadInfo  read1;
    SpanReadInfo read2;

    std::string readSeq1 = "GGTATTTTCGTCTGGGGGGT";
    std::string readSeq2 = "GGTATTTTCG";
    // when only mate-1 read available.
    bam_record bamRecord1;
    buildTestBamRecord(bamRecord1, 0, 200, 0, 270, 20, 15, "20M", readSeq1);
    bamRecord1.set_qname("Read-1");
    fragment1.read1.bamrec = bamRecord1;
    getFragInfo(fragment1, read1, read2);
    BOOST_REQUIRE(read1.isFwdStrand);
    BOOST_REQUIRE_EQUAL(read1.readSize, 20);
    BOOST_REQUIRE_EQUAL(read1.interval, GenomeInterval(0, 200, 220));
    // appoximate values of mate-2 as mate-2 is not available
    BOOST_REQUIRE(!read2.isFwdStrand);
    BOOST_REQUIRE_EQUAL(read2.readSize, 20);
    BOOST_REQUIRE_EQUAL(read2.interval, GenomeInterval(0, 270, 290));

    // When mate-1 and mate-2 both are available
    bam_record bamRecord2;
    buildTestBamRecord(bamRecord2, 0, 270, 0, 200, 20, 15, "10M", readSeq2);
    bamRecord2.toggle_is_fwd_strand();
    bamRecord2.toggle_is_mate_fwd_strand();
    bamRecord2.set_qname("Read-2");
    fragment1.read2.bamrec = bamRecord2;
    getFragInfo(fragment1, read1, read2);
    BOOST_REQUIRE(read1.isFwdStrand);
    BOOST_REQUIRE_EQUAL(read1.readSize, 20);
    BOOST_REQUIRE_EQUAL(read1.interval, GenomeInterval(0, 200, 220));
    BOOST_REQUIRE(!read2.isFwdStrand);
    BOOST_REQUIRE_EQUAL(read2.readSize, 10);
    BOOST_REQUIRE_EQUAL(read2.interval, GenomeInterval(0, 270, 280));

    // When only mate-2 read available
    SVCandidateSetSequenceFragment fragment2;
    fragment2.read2.bamrec = bamRecord2;
    getFragInfo(fragment2, read1, read2);
    BOOST_REQUIRE(!read2.isFwdStrand);
    BOOST_REQUIRE_EQUAL(read2.readSize, 10);
    BOOST_REQUIRE_EQUAL(read2.interval, GenomeInterval(0, 270, 280));
    // approximate values of mate-1 as mate-1 is not available.
    BOOST_REQUIRE(read1.isFwdStrand);
    BOOST_REQUIRE_EQUAL(read1.readSize, 10);
    BOOST_REQUIRE_EQUAL(read1.interval, GenomeInterval(0, 200, 210));
}

// Test the terminal information. Following cases need to tested
// 1. Whether read orientation information has copied properly or not
// 2. Whether read size information has copied properly or not
// 3. Whether chromosome information has copied properly or not
// 4. If orientation if forward then terminal position should be start
//    of the interval else terminal position is end of the interval
BOOST_AUTO_TEST_CASE( test_getTerminal )
{
    SpanReadInfo readInfo;
    SpanTerminal terminal;
    // case of forward strand
    readInfo.isFwdStrand = true;
    readInfo.readSize = 150;
    readInfo.interval = GenomeInterval(0, 100, 250);
    getTerminal(readInfo, terminal);
    BOOST_REQUIRE(terminal.isFwdStrand);
    BOOST_REQUIRE_EQUAL(terminal.readSize, 150);
    BOOST_REQUIRE_EQUAL(terminal.pos, 100);
    BOOST_REQUIRE_EQUAL(terminal.tid, 0);

    // case for reverse strand
    readInfo.isFwdStrand = false;
    getTerminal(readInfo, terminal);
    BOOST_REQUIRE_EQUAL(terminal.pos, 250);
}

// Test the read fragment probability if that fragment supports SV.
// FragProb = min(C, 1-C),
// C = fragDist.cdf(frag1Size + frag2Size), where
//            fragDist = sample fragment distribution,
//            cdf = cumulative density function
//            frag1Size = distance from center pos of BP1 to read1
//            frag2Size = distance from center pos of BP2 to read2
// Following cases need to be tested:
// 1. If a fragment supports SV frag probability is calculated as mentioned above
// 2. Fragment probability should be greater than 0.0001f.
// 3. Distance from breakpoint center position to read position should be greater than 50
// 4. For RNA FragProb = max(FragProb, 0.0001f)
// 5. If read pair is inter chromosomal but breakpoint pair is not inter chromosomal, it throws an exception.
// 6. If Read1 chromosome id and bp1 chromosome id is not same, it throws an exception
// 7. If Read2 chromosome id and bp2 chromosome id is not same, it throws an exception
// 8. if Read1 is in forward strand and sv.bp1.state != SVBreakendState::RIGHT_OPEN, it throws an exception
// 9. If Read2 is in reverse strand and sv.bp1.state == SVBreakendState::RIGHT_OPEN, it throws an exception
BOOST_AUTO_TEST_CASE( test_getFragProb )
{
    PairOptions options1(false); // false means it is a DNA sample
    bool isFragSupportSV;
    float fragProb;
    SVCandidate candidate;
    // bp1 center pos = 254
    candidate.bp1.interval = GenomeInterval(0, 250, 260);
    candidate.bp1.state = SVBreakendState::RIGHT_OPEN;
    // bp2 center pos = 276
    candidate.bp2.interval = GenomeInterval(0, 270, 280);
    candidate.bp2.state = SVBreakendState::LEFT_OPEN;
    SVCandidateSetSequenceFragment fragment;
    // creating fragment distribution
    SizeDistribution fragDistribution;
    for (unsigned i(0); i < 250; ++i)
    {
        fragDistribution.addObservation(50);
        fragDistribution.addObservation(75);
        fragDistribution.addObservation(100);
        fragDistribution.addObservation(125);
    }
    std::string queryseq1 = "AGCTGACTGATCGATTTTTTACGTAGAGGAGCTTTGACGTATGAGCCTGATATGAGCCTG";
    std::string queryseq2 = "TGACGTATGAGCCTGATATGAGCCT";

    // This is mate-1 read which is in forward strand
    bam_record bamRecord1;
    buildTestBamRecord(bamRecord1, 0, 200, 0, 300, 125, 15, "60M", queryseq1);
    bamRecord1.set_qname("Read-1");
    fragment.read1.bamrec = bamRecord1;

    // This is mate-2 read which is in reverse strand
    bam_record bamRecord2;
    buildTestBamRecord(bamRecord2, 0, 300, 0, 200, 125, 15, "25M", queryseq2);
    bamRecord2.toggle_is_fwd_strand();
    bamRecord2.toggle_is_mate_fwd_strand();
    bamRecord2.set_qname("Read-2");
    fragment.read2.bamrec = bamRecord2;


    float expected = 0.25;
    // This is the positive test case of DNA sample.
    // Only Case-1, Case-2 and case-3 are satisfied.
    // Orientation of reads and breakpoints are fine
    getFragProb(options1, candidate, fragment, fragDistribution, true, isFragSupportSV, fragProb);
    BOOST_REQUIRE(isFragSupportSV);
    BOOST_REQUIRE_EQUAL(((int)(fragProb * 100 + 0.5)) / 100.0, expected);

    // This is similar test case as above but here fragment probality is performed
    // on RNA sample. So case-4 is also considered here.
    PairOptions options2(true);
    getFragProb(options2, candidate, fragment, fragDistribution, true, isFragSupportSV, fragProb);
    BOOST_REQUIRE(isFragSupportSV);
    BOOST_REQUIRE_EQUAL(((int)(fragProb * 100 + 0.5)) / 100.0, expected);


    // BamRecord-1 and bamRecord-2 both are in forward strand
    // If both are in forward strand, then it will check which breakpoint is closer
    // to which read.
    bamRecord1.toggle_is_mate_fwd_strand();
    fragment.read1.bamrec = bamRecord1;
    bamRecord2.toggle_is_fwd_strand();
    bamRecord2.toggle_is_mate_fwd_strand();
    fragment.read2.bamrec = bamRecord2;
    // bp2 center pos = 224
    candidate.bp2.interval = GenomeInterval(0, 220, 230);
    // case-3 is not satisfied here.
    // position of bamRecord-1 is less than position of bamRecord-2 and also bp2 center position
    // less than bp1 center position. So distance calculation will happen between (bamRecord-1 and bp2)
    // and between (bamRecord-2 and bp1).
    // So Read-1 distance = 224 - 200 = 24 which is less than < 50
    getFragProb(options1, candidate, fragment, fragDistribution, true, isFragSupportSV, fragProb);
    BOOST_REQUIRE(!isFragSupportSV);
    BOOST_REQUIRE_EQUAL(fragProb, 0);


    // Interchromosomal read pair
    bam_record bamRecord3;
    buildTestBamRecord(bamRecord3, 0, 200, 1, 300, -1, 15, "60M", queryseq1);
    bamRecord3.set_qname("Read-3");
    fragment.read1.bamrec = bamRecord3;
    bam_record bamRecord4;
    buildTestBamRecord(bamRecord4, 1, 300, 0, 200, -1, 15, "25M", queryseq2);
    bamRecord4.toggle_is_fwd_strand();
    bamRecord4.toggle_is_mate_fwd_strand();
    fragment.read2.bamrec = bamRecord4;
    // Breakpoint pair is not interchromosomal
    candidate.bp2.interval = GenomeInterval(0, 270, 280);
    // Case-5 are designed here
    BOOST_CHECK_THROW(getFragProb(options1, candidate, fragment, fragDistribution, true, isFragSupportSV, fragProb), illumina::common::GeneralException);

    // Read pair and breakpoint pair are inter chromosomal, but
    // chromosome of read-1 and BP-1 are not matching. Case-6 is designed here.
    candidate.bp1.interval = GenomeInterval(2, 250, 260);
    candidate.bp2.interval = GenomeInterval(0, 270, 280);
    BOOST_CHECK_THROW(getFragProb(options1, candidate, fragment, fragDistribution, true, isFragSupportSV, fragProb), illumina::common::GeneralException);

    // Read pair and breakpoint pair are inter chromosomal, but
    // chromosome of read-2 and BP-2 are not matching. Case-7 is designed here.
    candidate.bp1.interval = GenomeInterval(0, 250, 260);
    candidate.bp2.interval = GenomeInterval(2, 270, 280);
    BOOST_CHECK_THROW(getFragProb(options1, candidate, fragment, fragDistribution, true, isFragSupportSV, fragProb), illumina::common::GeneralException);

    // Case-8 is designed here where read orientation and breakpoint orientation are not matching
    bamRecord3.toggle_is_fwd_strand();
    bamRecord3.get_data()->core.mtid = 0;
    fragment.read1.bamrec = bamRecord3;
    bamRecord4.set_target_id(0);
    fragment.read2.bamrec = bamRecord4;
    candidate.bp1.interval = GenomeInterval(0, 250, 260);
    candidate.bp1.state = SVBreakendState::RIGHT_OPEN;
    candidate.bp2.interval = GenomeInterval(0, 270, 280);
    candidate.bp2.state = SVBreakendState::RIGHT_OPEN;
    BOOST_CHECK_THROW(getFragProb(options1, candidate, fragment, fragDistribution, true, isFragSupportSV, fragProb), illumina::common::GeneralException);

    // Case-9 is designed here where read orientation and breakpoint orientation are not matching
    bamRecord3.toggle_is_fwd_strand();
    fragment.read1.bamrec = bamRecord3;
    BOOST_CHECK_THROW(getFragProb(options1, candidate, fragment, fragDistribution, true, isFragSupportSV, fragProb), illumina::common::GeneralException);
}

// Test whether a bam read satisfies the following criteria
// 1. Start of bam read should overlap with the search range
// 2. Fragment length should be greater than or equal to min fragment length of the sample
// 3. Fragment length should be less than or equal to max fragment length of the sample
// 4. minimum of (breakpoint center pos - fragment start + 1) and (fragment end - breakpoint center pos)
//    should be greater than minimum fragment threshold which is 50
BOOST_AUTO_TEST_CASE( test_processBamProcList )
{
    // whether alignment file tumor or normal
    const std::vector<bool> bamFileInfo = {false};
    const bam_header_info bamHeader(buildTestBamHeader());
    std::unique_ptr<SVLocusScanner> scanner(buildTestSVLocusScanner(bamHeader));
    const PairOptions options1(false);
    SVCandidate candidate;
    candidate.insertSeq = "GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGT"
                          "ATTTTCGTCTGGGGGGTGTGCACGCGATAGCATTGCGAGACGCTGGA";
    // BP1 center pos = 159
    candidate.bp1.interval.range = known_pos_range2(100, 220);
    // BP2 center pos = 309
    candidate.bp2.interval.range = known_pos_range2(250,370);
    SVEvidence evidence;
    evidence.samples.resize(1);
    std::shared_ptr<SVScorePairRefProcessor> processor(new SVScorePairRefProcessor(bamFileInfo,
                                                       scanner.operator*(), options1, candidate, true, evidence));
    std::vector<SVScorer::pairProcPtr> pairProcList = {processor};
    std::vector<SVScorer::streamPtr> bamStreams = {bamStream};
    SVId id;
    SupportSamples suppSamples;
    suppSamples.supportSamples.resize(1);
    processBamProcList(bamStreams, id, pairProcList, suppSamples);
    // bam read start = 9. It is not overlapping with search range [84, 235).
    // As a result of this, fragment support of this bam read will not be set.
    BOOST_REQUIRE(!evidence.getSampleEvidence(0)["bamRecord1"].ref.bp1.isFragmentSupport);
    // Bam read start = 109. It is overlapping with the search range. But
    // its fragment length is 49 which is less than minimum fragment length(50) of the sample
    // As a result of this, fragment support of this bam read will not be set.
    BOOST_REQUIRE(!evidence.getSampleEvidence(0)["bamRecord2"].ref.bp1.isFragmentSupport);
    // Bam read start = 109. It is overlapping with the search range. But
    // its fragment length is 130 which is greater than maximum fragment length(125) of the
    // sample. As a result of this, fragment support of this bam read will not be set.
    BOOST_REQUIRE(!evidence.getSampleEvidence(0)["bamRecord3"].ref.bp1.isFragmentSupport);
    // Here point-4 is not satisfied. As a result of this, fragment support of this
    // bam read will not be set.
    BOOST_REQUIRE(evidence.getSampleEvidence(0)["bamRecord4"].ref.bp1.isFragmentSupport);
}

// count the read pairs supporting the alternate allele in each sample.
// Also for each supporting fragment it will compute the fragment probability as explained
// in test_getFragProb.
BOOST_AUTO_TEST_CASE( test_processExistingAltPairInfo )
{
    const bam_header_info bamHeader(buildTestBamHeader());
    std::unique_ptr<SVLocusScanner> scanner(buildTestSVLocusScanner(bamHeader));
    GSCOptions options;
    options.alignFileOpt.alignmentFilenames = {bamFileName};
    SVScorer scorer(options, scanner.operator*(), bamHeader);
    TestSVScorer fSVScorer;
    PairOptions pairOptions(false); // false means it is a DNA sample

    // Candidate SV
    SVCandidate candidate;
    candidate.bp1.interval = GenomeInterval(0, 250, 260);
    candidate.bp1.state = SVBreakendState::RIGHT_OPEN;
    candidate.bp2.interval = GenomeInterval(0, 270, 280);
    candidate.bp2.state = SVBreakendState::LEFT_OPEN;

    std::string queryseq1 = "AGCTGACTGATCGATTTTTTACGTAGAGGAGCTTTGACGTATGAGCCTGATATGAGCCTG";
    std::string queryseq2 = "TGACGTATGAGCCTGATATGAGCCT";
    bam_record bamRecord1;
    buildTestBamRecord(bamRecord1, 0, 200, 0, 300, 125, 20, "60M", queryseq1);
    bamRecord1.set_qname("Read-1");
    bamRecord1.toggle_is_first();
    bam_record bamRecord2;
    buildTestBamRecord(bamRecord2, 0, 300, 0, 200, 125, 10, "25M", queryseq2);
    bamRecord2.toggle_is_second();
    bamRecord2.toggle_is_fwd_strand();
    bamRecord2.toggle_is_mate_fwd_strand();
    bamRecord2.set_qname("Read-1");
    float expected = 0.25;

    // Creating a fragment
    SVCandidateSetData candidateSetData;
    SVCandidateSetSequenceFragmentSampleGroup& group = candidateSetData.getDataGroup(0);
    group.add(bamHeader, bamRecord1, false, true, true);
    group.add(bamHeader, bamRecord2, false, true, true);
    SVSequenceFragmentAssociation association(0, SVEvidenceType::PAIR);
    group.begin().operator*().svLink.push_back(association);
    SVId id;
    SVEvidence evidence;
    SupportSamples suppSamples;
    suppSamples.supportSamples.resize(1);
    evidence.samples.resize(1);
    fSVScorer.processExistingAltPairInfo(scorer, pairOptions, candidateSetData, candidate, id, evidence, suppSamples);

    for (const SVCandidateSetSequenceFragment& fragment : group)
    {
        SVFragmentEvidence fragmentEvidence(evidence.getSampleEvidence(0)[fragment.qname()]);
        SVFragmentEvidenceAllele alt(fragmentEvidence.alt);
        // Whether read-1 is scanned
        BOOST_REQUIRE(fragmentEvidence.read1.isScanned);
        // Whether read-2 is scanned
        BOOST_REQUIRE(fragmentEvidence.read2.isScanned);
        // Read-1 is anchored read as mapping quality (20) greater than min mappig quality (15)
        BOOST_REQUIRE(fragmentEvidence.read1.isAnchored(false));
        // Read-1 is anchored read as mapping quality (10) less than min mappig quality(15)
        BOOST_REQUIRE(!fragmentEvidence.read2.isAnchored(false));
        // fragment pronbability as explained in test_getFragProb
        BOOST_REQUIRE_EQUAL(((int)(alt.bp1.fragLengthProb * 100 + 0.5)) / 100.0, expected);
        BOOST_REQUIRE_EQUAL(((int)(alt.bp2.fragLengthProb * 100 + 0.5)) / 100.0, expected);
        // Fragments are supporting this breakpoints as explained in test_processBamProcList
        BOOST_REQUIRE(alt.bp1.isFragmentSupport);
        BOOST_REQUIRE(alt.bp2.isFragmentSupport);
    }

}

// Following two cases are designed here:
// Let's say, fragment length in extreme 5th-95th percentiles over all read groups is F and
// deletion size from candidate SV is D.
// 1. If D <= 2*F then it is incomplete alt pair info, so it should not support the breakpoint
// 2. If D > 2*F, then fragment probability is cacluated as explained in test_getFragProb
BOOST_AUTO_TEST_CASE( test_getSVPairSupport )
{
    // Dummy Assembly data is created
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
    candidateAssemblyData.isCandidateSpanning = true;
    // Test SV  locus created with max fragment length in 5th to 95th percentile is 125
    const bam_header_info bamHeader(buildTestBamHeader());
    std::unique_ptr<SVLocusScanner> scanner(buildTestSVLocusScanner(bamHeader));

    GSCOptions options;
    options.alignFileOpt.alignmentFilenames = { bamFileName };
    SVScorer scorer(options, scanner.operator*(), bamHeader);
    TestSVScorer fSVScorer;
    PairOptions pairOptions(false);
    SVCandidate candidate1;
    candidate1.bp1.interval = GenomeInterval(0, 250, 260);
    candidate1.bp1.state = SVBreakendState::RIGHT_OPEN;
    candidate1.bp2.interval = GenomeInterval(0, 270, 280);
    candidate1.bp2.state = SVBreakendState::LEFT_OPEN;
    std::string queryseq1 = "AGCTGACTGATCGATTTTTTACGTAGAGGAGCTTTGACGTATGAGCCTGATATGAGCCTG";
    std::string queryseq2 = "TGACGTATGAGCCTGATATGAGCCT";
    bam_record bamRecord1;
    buildTestBamRecord(bamRecord1, 0, 200, 0, 300, 125, 15, "60M", queryseq1);
    bamRecord1.set_qname("Read-1");
    bamRecord1.toggle_is_first();
    bam_record bamRecord2;
    buildTestBamRecord(bamRecord2, 0, 300, 0, 200, 125, 15, "25M", queryseq2);
    bamRecord2.toggle_is_second();
    bamRecord2.toggle_is_fwd_strand();
    bamRecord2.toggle_is_mate_fwd_strand();
    bamRecord2.set_qname("Read-1");

    // Creating fragment
    SVCandidateSetData candidateSetData1;
    SVCandidateSetSequenceFragmentSampleGroup& group1 = candidateSetData1.getDataGroup(0);
    group1.add(bamHeader, bamRecord1, false, true, true);
    group1.add(bamHeader, bamRecord2, false, true, true);
    SVSequenceFragmentAssociation association(0, SVEvidenceType::PAIR);
    group1.begin().operator*().svLink.push_back(association);
    SVId id;
    SVEvidence evidence1;
    SupportSamples suppSamples1;
    suppSamples1.supportSamples.resize(1);
    evidence1.samples.resize(1);
    candidate1.setPrecise();
    candidate1.assemblyAlignIndex = 0;
    fSVScorer.getSVPairSupport(scorer, candidateSetData1, candidateAssemblyData, candidate1, id, evidence1, suppSamples1);

    // According to description, F = 125 and D = 270-250 = 20
    // So D < 2*F (20 < 250). Nothing will be calculated.
    for (const SVCandidateSetSequenceFragment& fragment : group1)
    {
        SVFragmentEvidenceAllele alt(evidence1.getSampleEvidence(0)[fragment.qname()].alt);
        BOOST_REQUIRE_EQUAL(((int)(alt.bp1.fragLengthProb * 100 + 0.5)) / 100.0, 0);
        BOOST_REQUIRE_EQUAL(((int)(alt.bp2.fragLengthProb * 100 + 0.5)) / 100.0, 0);
        BOOST_REQUIRE(!alt.bp1.isFragmentSupport);
        BOOST_REQUIRE(!alt.bp2.isFragmentSupport);
    }

    SVCandidate candidate2;
    candidate2.bp1.interval = GenomeInterval(0, 250, 260);
    candidate2.bp1.state = SVBreakendState::RIGHT_OPEN;
    candidate2.bp2.interval = GenomeInterval(0, 770, 780);
    candidate2.bp2.state = SVBreakendState::LEFT_OPEN;
    bam_record bamRecord3;
    buildTestBamRecord(bamRecord3, 0, 200, 0, 800, 625, 15, "60M", queryseq1);
    bamRecord3.set_qname("Read-2");
    bamRecord3.toggle_is_first();

    bam_record bamRecord4;
    buildTestBamRecord(bamRecord4, 0, 800, 0, 200, 625, 15, "25M", queryseq2);
    bamRecord4.toggle_is_second();
    bamRecord4.toggle_is_fwd_strand();
    bamRecord4.toggle_is_mate_fwd_strand();
    bamRecord4.set_qname("Read-2");
    float expected = 0.25;

    SVCandidateSetData candidateSetData2;
    SVCandidateSetSequenceFragmentSampleGroup& group2 = candidateSetData2.getDataGroup(0);
    group2.add(bamHeader, bamRecord3, false, true, true);
    group2.add(bamHeader, bamRecord4, false, true, true);
    group2.begin().operator*().svLink.push_back(association);
    SVEvidence evidence2;
    SupportSamples suppSamples2;
    suppSamples2.supportSamples.resize(1);
    evidence2.samples.resize(1);
    candidate2.setPrecise();
    candidate2.assemblyAlignIndex = 0;
    fSVScorer.getSVPairSupport(scorer, candidateSetData2, candidateAssemblyData, candidate2, id, evidence2, suppSamples2);
    // Here F = 125, But D = 770-250 = 520 which is greater than 2*F(250).
    // Fragment probability is calculated as explained in test_getFragProb.
    for (const SVCandidateSetSequenceFragment& fragment : group2)
    {
        SVFragmentEvidenceAllele alt(evidence2.getSampleEvidence(0)[fragment.qname()].alt);
        BOOST_REQUIRE_EQUAL(((int)(alt.bp1.fragLengthProb * 100 + 0.5)) / 100.0, expected);
        BOOST_REQUIRE_EQUAL(((int)(alt.bp2.fragLengthProb * 100 + 0.5)) / 100.0, expected);
        BOOST_REQUIRE(alt.bp1.isFragmentSupport);
        BOOST_REQUIRE(alt.bp2.isFragmentSupport);
    }
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()