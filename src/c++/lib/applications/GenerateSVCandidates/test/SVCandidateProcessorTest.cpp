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

#include <fstream>
#include <htsapi/vcf_streamer.hh>
#include "boost/test/unit_test.hpp"
#include "boost/make_unique.hpp"
#include "test/testAlignmentDataUtil.hh"
#include "test/testFileMakers.hh"
#include "manta/BamStreamerUtils.hh"
#include "test/testSVLocusScanner.hh"
#include "test/testSVLocusUtil.hh"

#include "options/SVLocusSetOptions.hh"
#include "svgraph/SVLocusSet.hh"
#include "SVCandidateProcessor.hh"
#include "SVCandidateProcessor.cpp"

/// The whole purpose of this test file is given an SV candidate, whether it is correctly writing
/// SV information in the corresponding vcf file for tumor, rna, diploid and somatic cases.

/// TestSVCandidateProcessor is a friend of GSCEdgeStatsManager. So that can access private
/// members of SVCandidateProcessor.
struct TestSVCandidateProcessor
{
    // closing all the streams, as we need to read those files before
    // program exits.
    void flushStreams(SVCandidateProcessor &candidateProcessor)
    {
        candidateProcessor._svWriter.tumfs.getStream().flush();
        candidateProcessor._svWriter.rnafs.getStream().flush();
        candidateProcessor._svWriter.dipfs.getStream().flush();
        candidateProcessor._svWriter.somfs.getStream().flush();
        candidateProcessor._svWriter.candfs.getStream().flush();
    }
};

// As for somatic, minimum two bam files(normal and tumor) are needed,
// so creating stats for two bam file.
std::unique_ptr<SVLocusScanner> buildSomaticSVLocusScanner(bam_header_info bamHeaderInfo)
{
    TestFilenameMaker filenameMaker;
    std:: string tmpFileName = filenameMaker.getFilename();
    ReadGroupStatsSet rstats;
    for (unsigned j(0); j < 2; j++)
    {
        ReadGroupLabel rgKey(("Bamfile" + std::to_string(j)).c_str(), std::to_string(j).c_str());
        ReadGroupStats rgStats;
        for (unsigned i(0); i < 250; ++i) {
            rgStats.fragStats.addObservation(50);
            rgStats.fragStats.addObservation(75);
            rgStats.fragStats.addObservation(100);
            rgStats.fragStats.addObservation(125);
        }
        rstats.setStats(rgKey, rgStats);
    }
    rstats.save(tmpFileName.c_str());
    ReadScannerOptions opts = ReadScannerOptions();
    opts.minCandidateVariantSize = 8;
    const ReadScannerOptions& constRefOpts(opts);
    TestAlignHeaderFileMaker alignFile(bamHeaderInfo);
    TestAlignHeaderFileMaker alignFile1(bamHeaderInfo);
    const std::vector<std::string> alignFilenameVector = { alignFile.getFilename(), alignFile1.getFilename()};
    return boost::make_unique<SVLocusScanner>(constRefOpts, tmpFileName,
                                              alignFilenameVector, false);
}

BOOST_AUTO_TEST_SUITE( SVCandidateProcessor_test_suite )

// Create Temporary bam streams of a bam file which contains
// three interchromosomal read pairs.
struct BamStream
{
    BamStream()
    {
        const bam_header_info bamHeader(buildTestBamHeader());

        std::string querySeq1 = "GTCTATCACCCTATTAACCACTCACGGGAGAAAAA";
        std::string querySeq2 = "AAAAATGCTCATCAGTTGATGATACGCCCGAGCAGATGCCAACACAGCAGCCATTCAAGAGTCA";
        bam_record bamRecord1;
        buildTestBamRecord(bamRecord1, 0, 8, 1, 69, -1, 15, "35M", querySeq1);
        bamRecord1.set_qname("Read-1");
        bamRecord1.toggle_is_first();
        bam_record bamRecord2;
        buildTestBamRecord(bamRecord2, 0, 8, 1, 69, -1, 15, "35M", querySeq1);
        bamRecord2.set_qname("Read-2");
        bamRecord2.toggle_is_first();
        bam_record bamRecord3;
        buildTestBamRecord(bamRecord3, 0, 8, 1, 69, -1, 15, "35M", querySeq1);
        bamRecord3.set_qname("Read-3");
        bamRecord3.toggle_is_first();
        bam_record bamRecord4;
        buildTestBamRecord(bamRecord4, 1, 69, 0, 8, -1, 50, "64M", querySeq2);
        bamRecord4.toggle_is_mate_fwd_strand();
        bamRecord4.toggle_is_fwd_strand();
        bamRecord4.set_qname("Read-1");
        bamRecord4.toggle_is_second();
        bam_record bamRecord5;
        buildTestBamRecord(bamRecord5, 1, 69, 0, 8, -1, 50, "64M", querySeq2);
        bamRecord5.set_qname("Read-2");
        bamRecord5.toggle_is_mate_fwd_strand();
        bamRecord5.toggle_is_fwd_strand();
        bamRecord5.toggle_is_second();
        bam_record bamRecord6;
        buildTestBamRecord(bamRecord6, 1, 69, 0, 8, -1, 50, "64M", querySeq2);
        bamRecord6.toggle_is_mate_fwd_strand();
        bamRecord6.toggle_is_fwd_strand();
        bamRecord6.set_qname("Read-3");
        bamRecord6.toggle_is_second();
        readsToAdd.push_back(bamRecord1);
        readsToAdd.push_back(bamRecord2);
        readsToAdd.push_back(bamRecord3);
        readsToAdd.push_back(bamRecord4);
        readsToAdd.push_back(bamRecord5);
        readsToAdd.push_back(bamRecord6);
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

BOOST_FIXTURE_TEST_SUITE( SVCandidateProcessor_test_suite, BamStream )

// A candidate sv needs to be checked after assembly whether it is good candidate sv or not
// A candidate sv is not good if either of the following three cases is  satisfied:
// 1. Post assembly spanning count is less than min candidate spanning count threshold (default 3).
// 2. Non-spanning low-res candidate went into assembly but did not produce a successful contig alignment
//    that means SV is imprecise.
// 3. Variant size is smaller than minCandidateVariantSize (default is 10)
// Successful cases are designed in test_tumorOnly, test_RNA, test_Diploid and test_Somatic.
BOOST_AUTO_TEST_CASE( test_JunctionFilter )
{
    const bam_header_info bamHeader(buildTestBamHeader());
    TestFilenameMaker fileMakerBase1;
    TestFilenameMaker fileMakerBase2;
    TestFilenameMaker fileMakerBase3;
    TestFilenameMaker fileMakerBase4;
    GSCOptions options;

    // Creating SV candidate
    SVCandidate candidate1;
    candidate1.setPrecise();
    candidate1.bp1.state = SVBreakendState::RIGHT_OPEN;
    candidate1.bp1.interval = GenomeInterval(0 , 40, 50);
    candidate1.bp1.lowresEvidence.add(0, 1);
    candidate1.bp2.state = SVBreakendState::LEFT_OPEN;
    candidate1.bp2.interval = GenomeInterval(1 , 65, 75);
    candidate1.bp2.lowresEvidence.add(0, 1);
    std::unique_ptr<SVLocusScanner> scanner(buildTestSVLocusScanner(bamHeader));
    std::string programName = "Manta";
    std::string version = "Test:0:1";
    EdgeInfo edgeInfo;
    edgeInfo.nodeIndex1 = 0;
    edgeInfo.nodeIndex2 = 1;
    SVMultiJunctionCandidate junctionCandidate;
    junctionCandidate.junction = {candidate1};
    // Adding fragment information
    SVCandidateSetData svData;
    SVCandidateSetSequenceFragmentSampleGroup& group = svData.getDataGroup(0);
    group.add(bamHeader, readsToAdd[0], false, true, true);
    group.add(bamHeader, readsToAdd[3], false, true, true);
    SVSequenceFragmentAssociation association(0, SVEvidenceType::PAIR);
    group.begin().operator*().svLink.push_back(association);
    SupportSamples svSupports;
    svSupports.supportSamples.resize(1);
    // Generate SV locus graph file
    // Two nodes are created.
    SVLocus locus1;
    locusAddPair(locus1, 0, 40, 50, 1, 65, 75);
    SVLocus locus2;
    locusAddPair(locus2, 1, 65, 75, 0, 40, 50);
    SVLocusSetOptions sopt;
    sopt.minMergeEdgeObservations = 1;
    SVLocusSet set1(sopt, bamHeader, {bamFileName});
    set1.merge(locus1);
    set1.merge(locus2);
    set1.checkState(true,true);
    TestFilenameMaker testFilenameMaker;
    std::string graphFilename(testFilenameMaker.getFilename());
    const char* testFilenamePtr(graphFilename.c_str());
    // serialize
    set1.save(testFilenamePtr);
    std::vector<bool> isInputJunctionFiltered;
    isInputJunctionFiltered.resize(1);
    SVCandidateAssemblyData candidateAssemblyData;
    candidateAssemblyData.isCandidateSpanning = true;
    std::vector <SVCandidateAssemblyData> assemblyData;
    assemblyData.push_back(candidateAssemblyData);
    // Case-1 is designed
    SVWriter writer(options, scanner.operator*(), bamHeader, programName.c_str(), version.c_str());
    writer.writeSV(edgeInfo, svData, assemblyData, junctionCandidate, isInputJunctionFiltered, svSupports);
    writer.candfs.getStream().flush();
    writer.tumfs.getStream().flush();
    std::ifstream candidateFile(options.candidateOutputFilename);
    std::string line;
    int count(0);

    // Post assembly spanning count is less than min candidate spanning count 3.
    if (candidateFile.is_open())
    {
        while (std::getline(candidateFile, line))
        {
            if (line.find("#") == std::string::npos)
            {
                count ++;
            }

        }
        candidateFile.close();
    }
    BOOST_REQUIRE_EQUAL(count, 0);

    // Case-2 is designed
    // in this case a non-spanning low-res candidate went into assembly but
    // did not produce a successful contig alignment:
    options.minCandidateSpanningCount = 1;
    SVCandidate candidate2;
    candidate2.bp1.state = SVBreakendState::COMPLEX;
    candidate2.bp1.interval = GenomeInterval(0 , 40, 50);
    candidate2.bp1.lowresEvidence.add(0, 1);
    candidate2.bp2.state = SVBreakendState::UNKNOWN;
    candidate2.bp2.interval = GenomeInterval(0 , 65, 75);
    candidate2.bp2.lowresEvidence.add(0, 1);
    count = 0;
    junctionCandidate.junction.clear();
    junctionCandidate.junction.push_back(candidate2);
    candidateAssemblyData.isCandidateSpanning = false;
    assemblyData.clear();
    assemblyData.push_back(candidateAssemblyData);
    writer.writeSV(edgeInfo, svData, assemblyData, junctionCandidate, isInputJunctionFiltered, svSupports);
    writer.candfs.getStream().flush();
    writer.tumfs.getStream().flush();
    std::ifstream candidateFile1(options.candidateOutputFilename);
    if (candidateFile1.is_open())
    {
        while (std::getline(candidateFile1, line))
        {
            if (line.find("#") == std::string::npos)
            {
                count ++;
            }

        }
        candidateFile1.close();
    }
    BOOST_REQUIRE_EQUAL(count, 0);

    // Case-3 is designed.
    // variant size is less than min variant size(40)
    options.minCandidateSpanningCount = 1;
    SVCandidate candidate3;
    candidate3.bp1.state = SVBreakendState::RIGHT_OPEN;
    candidate3.bp1.interval = GenomeInterval(0 , 40, 50);
    candidate3.bp1.lowresEvidence.add(0, 1);
    candidate3.bp2.state = SVBreakendState::LEFT_OPEN;
    candidate3.bp2.interval = GenomeInterval(0 , 65, 75);
    candidate3.bp2.lowresEvidence.add(0, 1);
    count = 0;
    junctionCandidate.junction.clear();
    junctionCandidate.junction.push_back(candidate2);
    candidateAssemblyData.isCandidateSpanning = true;
    assemblyData.clear();
    assemblyData.push_back(candidateAssemblyData);
    options.scanOpt.minCandidateVariantSize = 40;
    SVWriter writer2(options, scanner.operator*(), bamHeader, programName.c_str(), version.c_str());
    writer2.writeSV(edgeInfo, svData, assemblyData, junctionCandidate, isInputJunctionFiltered, svSupports);
    writer2.candfs.getStream().flush();
    writer2.tumfs.getStream().flush();
    std::ifstream candidateFile2(options.candidateOutputFilename);
    if (candidateFile2.is_open())
    {
        while (std::getline(candidateFile2, line))
        {
            if (line.find("#") == std::string::npos)
            {
                count ++;
            }

        }
        candidateFile1.close();
    }
    BOOST_REQUIRE_EQUAL(count, 0);
}

// Test the VCF records for tumor only caller
BOOST_AUTO_TEST_CASE( test_tumorOnly )
{
    const bam_header_info bamHeader(buildTestBamHeader());
    TestFilenameMaker fileMakerBase1;
    TestFilenameMaker fileMakerBase2;
    TestFilenameMaker fileMakerBase3;
    TestFilenameMaker fileMakerBase4;
    GSCOptions options;
    // Assembly options
    options.refineOpt.spanningAssembleOpt.minWordLength = 3;
    options.refineOpt.spanningAssembleOpt.maxWordLength = 9;
    options.refineOpt.spanningAssembleOpt.wordStepSize = 3;
    options.refineOpt.smallSVAssembleOpt.minWordLength = 3;
    options.refineOpt.smallSVAssembleOpt.maxWordLength = 9;
    options.refineOpt.smallSVAssembleOpt.wordStepSize = 3;
    // input/output files
    options.alignFileOpt.alignmentFilenames = {bamFileName};
    options.alignFileOpt.isAlignmentTumor = {true}; // only tumor bam file
    options.referenceFilename = getTestReferenceFilename();
    options.edgeRuntimeFilename = fileMakerBase1.getFilename();
    options.edgeStatsFilename = fileMakerBase2.getFilename();
    options.tumorOutputFilename = fileMakerBase3.getFilename();
    options.candidateOutputFilename = fileMakerBase4.getFilename();
    options.minCandidateSpanningCount = 1;
    TestStatsFileMaker statsFileMaker;
    options.statsFilename = statsFileMaker.getFilename();
    EdgeRuntimeTracker edgeTracker(options.edgeRuntimeFilename);
    GSCEdgeStatsManager edgeStatMan(options.edgeStatsFilename);

    // Creating SV candidate
    SVCandidate candidate1;
    candidate1.setPrecise();
    candidate1.bp1.state = SVBreakendState::RIGHT_OPEN;
    candidate1.bp1.interval = GenomeInterval(0 , 40, 50);
    candidate1.bp1.lowresEvidence.add(0, 1);
    candidate1.bp2.state = SVBreakendState::LEFT_OPEN;
    candidate1.bp2.interval = GenomeInterval(1 , 65, 75);
    candidate1.bp2.lowresEvidence.add(0, 1);
    SVCandidateAssemblyData candidateAssemblyData1;
    std::unique_ptr<SVLocusScanner> scanner(buildTestSVLocusScanner(bamHeader));
    std::string programName = "Manta";
    std::string version = "Test:0:1";
    EdgeInfo edgeInfo;
    edgeInfo.nodeIndex1 = 0;
    edgeInfo.nodeIndex2 = 1;
    SVMultiJunctionCandidate junctionCandidate;
    junctionCandidate.junction = {candidate1};
    std::vector <SVMultiJunctionCandidate> mjSvs;
    mjSvs.push_back(junctionCandidate);
    // Adding fragment information
    SVCandidateSetData svData;
    SVCandidateSetSequenceFragmentSampleGroup& group = svData.getDataGroup(0);
    group.add(bamHeader, readsToAdd[0], false, true, true);
    group.add(bamHeader, readsToAdd[3], false, true, true);
    SVSequenceFragmentAssociation association(0, SVEvidenceType::PAIR);
    group.begin().operator*().svLink.push_back(association);
    SupportSamples svSupports;
    svSupports.supportSamples.resize(1);
    // Generate SV locus graph file
    // Two nodes are created.
    SVLocus locus1;
    locusAddPair(locus1, 0, 40, 50, 1, 65, 75);
    SVLocus locus2;
    locusAddPair(locus2, 1, 65, 75, 0, 40, 50);
    SVLocusSetOptions sopt;
    sopt.minMergeEdgeObservations = 1;
    SVLocusSet set1(sopt, bamHeader, {bamFileName});
    set1.merge(locus1);
    set1.merge(locus2);
    set1.checkState(true,true);
    TestFilenameMaker testFilenameMaker;
    std::string graphFilename(testFilenameMaker.getFilename());
    const char* testFilenamePtr(graphFilename.c_str());
    // serialize
    set1.save(testFilenamePtr);
    SVCandidateProcessor candidateProcessor(options, scanner.operator*(), programName.c_str(),
                                            version.c_str(), set1, edgeTracker, edgeStatMan);
    candidateProcessor.evaluateCandidates(edgeInfo, mjSvs, svData, svSupports);
    TestSVCandidateProcessor testSVCandidateProcessor;
    testSVCandidateProcessor.flushStreams(candidateProcessor);
    std::ifstream tumorFile(options.tumorOutputFilename);
    std::string line;
    int count(0);
    // Excpected following tow vcf records:
    // Record-1:    chrFoo	41	MantaBND:0:0:1:0:0:0:0	C	C[chrBar:66[	.	.
    //              SVTYPE=BND;MATEID=MantaBND:0:0:1:0:0:0:1;CIPOS=0,9;HOMLEN=9;HOMSEQ=TCCATGCAT;BND_PAIR_COUNT=0;PAIR_COUNT=1
    // Record-2:    chrBar	66	MantaBND:0:0:1:0:0:0:1	T	]chrFoo:41]T	.	.
    //              SVTYPE=BND;MATEID=MantaBND:0:0:1:0:0:0:0;CIPOS=0,9;HOMLEN=9;HOMSEQ=ATCCTTCTT;BND_PAIR_COUNT=0;PAIR_COUNT=1
    if (tumorFile.is_open())
    {
        while (std::getline(tumorFile, line))
        {
            // ignoring all the headers
            if (line.find("#") == std::string::npos)
            {
                count ++;
            }
            if (count == 1)
                BOOST_REQUIRE(line.find("C[chrBar:66[") != std::string::npos);
            else if (count == 2)
                BOOST_REQUIRE(line.find("]chrFoo:41]T") != std::string::npos);

        }
        tumorFile.close();
    }
    BOOST_REQUIRE_EQUAL(count, 2);
}

// Test the VCF records for RNA caller
BOOST_AUTO_TEST_CASE( test_RNA )
{
    const bam_header_info bamHeader(buildTestBamHeader());
    TestFilenameMaker fileMakerBase1;
    TestFilenameMaker fileMakerBase2;
    TestFilenameMaker fileMakerBase3;
    TestFilenameMaker filenameMaker4;
    GSCOptions options;
    // Assembly options
    options.refineOpt.spanningAssembleOpt.minWordLength = 3;
    options.refineOpt.spanningAssembleOpt.maxWordLength = 9;
    options.refineOpt.spanningAssembleOpt.wordStepSize = 3;
    options.refineOpt.smallSVAssembleOpt.minWordLength = 3;
    options.refineOpt.smallSVAssembleOpt.maxWordLength = 9;
    options.refineOpt.smallSVAssembleOpt.wordStepSize = 3;
    // Input/output
    options.alignFileOpt.alignmentFilenames = {bamFileName};
    options.alignFileOpt.isAlignmentTumor = {false}; // Not tumor
    options.isRNA = true;// This is an rna sample
    options.referenceFilename = getTestReferenceFilename();
    options.edgeRuntimeFilename = fileMakerBase1.getFilename();
    options.edgeStatsFilename = fileMakerBase2.getFilename();
    options.rnaOutputFilename = fileMakerBase3.getFilename();
    options.candidateOutputFilename = filenameMaker4.getFilename();
    options.minCandidateSpanningCount = 1;
    options.minScoredVariantSize = 20;
    TestStatsFileMaker statsFileMaker;
    options.statsFilename = statsFileMaker.getFilename();
    EdgeRuntimeTracker edgeTracker(options.edgeRuntimeFilename);
    GSCEdgeStatsManager edgeStatMan(options.edgeStatsFilename);
    // Creating SV Candidate
    SVCandidate candidate1;
    candidate1.setPrecise();
    candidate1.bp1.state = SVBreakendState::RIGHT_OPEN;
    candidate1.bp1.interval = GenomeInterval(0 , 40, 50);
    candidate1.bp1.lowresEvidence.add(0, 1);
    candidate1.bp2.state = SVBreakendState::LEFT_OPEN;
    candidate1.bp2.interval = GenomeInterval(1 , 65, 75);
    candidate1.bp2.lowresEvidence.add(0, 1);
    SVCandidateAssemblyData candidateAssemblyData1;
    std::unique_ptr<SVLocusScanner> scanner(buildTestSVLocusScanner(bamHeader));
    std::string programName = "Manta";
    std::string version = "Test:0:1";
    EdgeInfo edgeInfo;
    edgeInfo.nodeIndex1 = 0;
    edgeInfo.nodeIndex2 = 1;
    SVMultiJunctionCandidate junctionCandidate;
    junctionCandidate.junction = {candidate1};
    std::vector <SVMultiJunctionCandidate> mjSvs;
    mjSvs.push_back(junctionCandidate);
    // Adding fragment data
    SVCandidateSetData svData;
    SVCandidateSetSequenceFragmentSampleGroup& group = svData.getDataGroup(0);
    group.add(bamHeader, readsToAdd[0], false, true, true);
    group.add(bamHeader, readsToAdd[3], false, true, true);
    SVSequenceFragmentAssociation association(0, SVEvidenceType::PAIR);
    group.begin().operator*().svLink.push_back(association);
    SupportSamples svSupports;
    svSupports.supportSamples.resize(1);
    // Generate SV locus graph file
    // Two nodes are created.
    SVLocus locus1;
    locusAddPair(locus1, 0, 40, 50, 1, 65, 75);
    SVLocus locus2;
    locusAddPair(locus2, 1, 65, 75, 0, 40, 50);
    SVLocusSetOptions sopt;
    sopt.minMergeEdgeObservations = 10;
    SVLocusSet set1(sopt, bamHeader, {bamFileName});
    set1.merge(locus1);
    set1.merge(locus2);
    set1.checkState(true,true);
    TestFilenameMaker testFilenameMaker;
    std::string graphFilename = testFilenameMaker.getFilename();
    const char* testFilenamePtr(graphFilename.c_str());
    // serialize
    set1.save(testFilenamePtr);
    SVCandidateProcessor candidateProcessor(options, scanner.operator*(), programName.c_str(),
                                            version.c_str(), set1, edgeTracker, edgeStatMan);
    candidateProcessor.evaluateCandidates(edgeInfo, mjSvs, svData, svSupports);
    TestSVCandidateProcessor testSVCandidateProcessor;
    testSVCandidateProcessor.flushStreams(candidateProcessor);

    std::ifstream rnaFile(options.rnaOutputFilename);
    std::string line;
    int count(0);
    // Following records are expected:
    // Record-1:    chrFoo	41	MantaBND:0:0:1:0:0:0:0	C	C[chrBar:66[	.	.
    //              SVTYPE=BND;MATEID=MantaBND:0:0:1:0:0:0:1;CIPOS=0,9;HOMLEN=9;HOMSEQ=TCCATGCAT;BND_PAIR_COUNT=0;PAIR_COUNT=1
    // Record-2:    chrBar	66	MantaBND:0:0:1:0:0:0:1	T	]chrFoo:41]T	.	.
    //              SVTYPE=BND;MATEID=MantaBND:0:0:1:0:0:0:0;CIPOS=0,9;HOMLEN=9;HOMSEQ=ATCCTTCTT;BND_PAIR_COUNT=0;PAIR_COUNT=1
    if (rnaFile.is_open())
    {
        while (std::getline(rnaFile, line))
        {
            if (line.find("#") == std::string::npos)
            {
                count ++;
            }
            if (count == 1)
                BOOST_REQUIRE(line.find("C[chrBar:66[") != std::string::npos);
            else if (count == 2)
                BOOST_REQUIRE(line.find("]chrFoo:41]T") != std::string::npos);

        }
        rnaFile.close();
    }
    BOOST_REQUIRE_EQUAL(count, 2);
}

// Test the VCF records for Diploid caller
BOOST_AUTO_TEST_CASE( test_Diploid )
{
    const bam_header_info bamHeader(buildTestBamHeader());
    TestFilenameMaker fileMakerBase1;
    TestFilenameMaker fileMakerBase2;
    TestFilenameMaker fileMakerBase3;
    TestFilenameMaker filenameMaker4;
    GSCOptions options;

    // Assembly options
    options.refineOpt.spanningAssembleOpt.minWordLength = 3;
    options.refineOpt.spanningAssembleOpt.maxWordLength = 9;
    options.refineOpt.spanningAssembleOpt.wordStepSize = 3;
    options.refineOpt.smallSVAssembleOpt.minWordLength = 3;
    options.refineOpt.smallSVAssembleOpt.maxWordLength = 9;
    options.refineOpt.smallSVAssembleOpt.wordStepSize = 3;
    options.alignFileOpt.alignmentFilenames = {bamFileName};
    options.alignFileOpt.isAlignmentTumor = {false};
    options.referenceFilename = getTestReferenceFilename();
    options.edgeRuntimeFilename = fileMakerBase1.getFilename();
    options.edgeStatsFilename = fileMakerBase2.getFilename();
    options.diploidOutputFilename = fileMakerBase3.getFilename();
    options.candidateOutputFilename = filenameMaker4.getFilename();
    options.minCandidateSpanningCount = 1;
    options.minScoredVariantSize = 20;
    options.diploidOpt.minOutputAltScore = 0;

    TestStatsFileMaker statsFileMaker;
    options.statsFilename = statsFileMaker.getFilename();
    EdgeRuntimeTracker edgeTracker(options.edgeRuntimeFilename);
    GSCEdgeStatsManager edgeStatMan(options.edgeStatsFilename);

    // Creating SV candidate
    SVCandidate candidate1;
    candidate1.setPrecise();
    candidate1.bp1.state = SVBreakendState::RIGHT_OPEN;
    candidate1.bp1.interval = GenomeInterval(0 , 40, 50);
    candidate1.bp1.lowresEvidence.add(0, 1);
    candidate1.bp2.state = SVBreakendState::LEFT_OPEN;
    candidate1.bp2.interval = GenomeInterval(1 , 65, 75);
    candidate1.bp2.lowresEvidence.add(0, 1);
    SVCandidateAssemblyData candidateAssemblyData1;
    std::unique_ptr<SVLocusScanner> scanner(buildTestSVLocusScanner(bamHeader));
    std::string programName = "Manta";
    std::string version = "Test:0:1";
    EdgeInfo edgeInfo;
    edgeInfo.nodeIndex1 = 0;
    edgeInfo.nodeIndex2 = 1;
    SVMultiJunctionCandidate junctionCandidate;
    junctionCandidate.junction = {candidate1};
    std::vector <SVMultiJunctionCandidate> mjSvs;
    mjSvs.push_back(junctionCandidate);
    // Adding fragment information
    SVCandidateSetData svData;
    SVCandidateSetSequenceFragmentSampleGroup& group = svData.getDataGroup(0);
    group.add(bamHeader, readsToAdd[0], false, true, true);
    group.add(bamHeader, readsToAdd[3], false, true, true);
    SVSequenceFragmentAssociation association(0, SVEvidenceType::PAIR);
    group.begin().operator*().svLink.push_back(association);
    SupportSamples svSupports;
    svSupports.supportSamples.resize(1);
    // Generate SV locus graph file
    // Two nodes are created.
    SVLocus locus1;
    locusAddPair(locus1, 0, 40, 50, 1, 65, 75);
    SVLocus locus2;
    locusAddPair(locus2, 1, 65, 75, 0, 40, 50);
    SVLocusSetOptions sopt;
    sopt.minMergeEdgeObservations = 10;
    SVLocusSet set1(sopt, bamHeader, {bamFileName});
    set1.merge(locus1);
    set1.merge(locus2);
    set1.checkState(true,true);
    TestFilenameMaker testFilenameMaker;
    std::string graphFilename = testFilenameMaker.getFilename();
    const char* testFilenamePtr(graphFilename.c_str());
    // serialize
    set1.save(testFilenamePtr);
    SVCandidateProcessor candidateProcessor(options, scanner.operator*(), programName.c_str(),
                                            version.c_str(), set1, edgeTracker, edgeStatMan);
    candidateProcessor.evaluateCandidates(edgeInfo, mjSvs, svData, svSupports);
    TestSVCandidateProcessor testSVCandidateProcessor;
    testSVCandidateProcessor.flushStreams(candidateProcessor);

    std::ifstream diploidFile(options.diploidOutputFilename);
    std::string line;
    int count(0);
    // Record-1:    chrFoo	41	MantaBND:0:0:1:0:0:0:0	C	C[chrBar:66[	.	.
    //              SVTYPE=BND;MATEID=MantaBND:0:0:1:0:0:0:1;CIPOS=0,9;HOMLEN=9;HOMSEQ=TCCATGCAT;BND_PAIR_COUNT=0;PAIR_COUNT=1
    // Record-2:    chrBar	66	MantaBND:0:0:1:0:0:0:1	T	]chrFoo:41]T	.	.
    //              SVTYPE=BND;MATEID=MantaBND:0:0:1:0:0:0:0;CIPOS=0,9;HOMLEN=9;HOMSEQ=ATCCTTCTT;BND_PAIR_COUNT=0;PAIR_COUNT=1
    if (diploidFile.is_open())
    {
        while (std::getline(diploidFile, line))
        {
            if (line.find("#") == std::string::npos)
            {
                count ++;
            }
            if (count == 1)
                BOOST_REQUIRE(line.find("C[chrBar:66[") != std::string::npos);
            else if (count == 2)
                BOOST_REQUIRE(line.find("]chrFoo:41]T") != std::string::npos);

        }
        diploidFile.close();
    }
    BOOST_REQUIRE_EQUAL(count, 2);
}

// Test the VCF records for somatic caller
BOOST_AUTO_TEST_CASE( test_Somatic )
{
    const bam_header_info bamHeader(buildTestBamHeader());
    TestFilenameMaker fileMakerBase1;
    TestFilenameMaker fileMakerBase2;
    TestFilenameMaker fileMakerBase3;
    TestFilenameMaker fileMakerBase4;
    TestFilenameMaker filenameMaker5;
    GSCOptions options;

    // Assembly options
    options.refineOpt.spanningAssembleOpt.minWordLength = 3;
    options.refineOpt.spanningAssembleOpt.maxWordLength = 9;
    options.refineOpt.spanningAssembleOpt.wordStepSize = 3;
    options.refineOpt.smallSVAssembleOpt.minWordLength = 3;
    options.refineOpt.smallSVAssembleOpt.maxWordLength = 9;
    options.refineOpt.smallSVAssembleOpt.wordStepSize = 3;

    // Input/output
    options.alignFileOpt.alignmentFilenames = {bamFileName, bamFileName};
    options.alignFileOpt.isAlignmentTumor = {false, true}; // normal and tumor
    options.referenceFilename = getTestReferenceFilename();
    options.edgeRuntimeFilename = fileMakerBase1.getFilename();
    options.edgeStatsFilename = fileMakerBase2.getFilename();
    options.diploidOutputFilename = fileMakerBase3.getFilename();
    options.somaticOutputFilename = fileMakerBase4.getFilename();
    options.candidateOutputFilename = filenameMaker5.getFilename();

    TestFilenameMaker depthFileNameMaker;
    buildTestChromosomeDepthFile(depthFileNameMaker.getFilename());
    options.chromDepthFilename = depthFileNameMaker.getFilename();
    options.minCandidateSpanningCount = 1;
    options.minScoredVariantSize = 20;
    options.diploidOpt.minOutputAltScore = 0;
    options.somaticOpt.minOutputSomaticScore = 0;
    TestStatsFileMaker statsFileMaker;
    options.statsFilename = statsFileMaker.getFilename();
    EdgeRuntimeTracker edgeTracker(options.edgeRuntimeFilename);
    GSCEdgeStatsManager edgeStatMan(options.edgeStatsFilename);
    // Creating SV candidates
    SVCandidate candidate1;
    candidate1.setPrecise();
    candidate1.bp1.state = SVBreakendState::RIGHT_OPEN;
    candidate1.bp1.interval = GenomeInterval(0 , 40, 50);
    candidate1.bp1.lowresEvidence.add(0, 9);
    candidate1.bp1.lowresEvidence.add(1, 9);
    candidate1.bp2.state = SVBreakendState::LEFT_OPEN;
    candidate1.bp2.interval = GenomeInterval(1 , 65, 75);
    candidate1.bp2.lowresEvidence.add(0, 9);
    candidate1.bp2.lowresEvidence.add(1, 9);

    SVCandidate candidate2;
    candidate2.setPrecise();
    candidate2.bp1.state = SVBreakendState::RIGHT_OPEN;
    candidate2.bp1.interval = GenomeInterval(0 , 10, 20);
    candidate2.bp1.lowresEvidence.add(0, 1);
    candidate2.bp2.state = SVBreakendState::RIGHT_OPEN;
    candidate2.bp2.interval = GenomeInterval(1 , 100, 110);
    candidate2.bp2.lowresEvidence.add(0, 1);

    SVCandidateAssemblyData candidateAssemblyData1;
    std::unique_ptr<SVLocusScanner> scanner(buildSomaticSVLocusScanner(bamHeader));
    std::string programName = "Manta";
    std::string version = "123";
    EdgeInfo edgeInfo;
    edgeInfo.nodeIndex1 = 0;
    edgeInfo.nodeIndex2 = 1;
    SVMultiJunctionCandidate junctionCandidate;
    junctionCandidate.junction = {candidate1, candidate2};
    std::vector <SVMultiJunctionCandidate> mjSvs;
    mjSvs.push_back(junctionCandidate);
    // Adding fragment data;
    SVCandidateSetData svData;
    SVCandidateSetSequenceFragmentSampleGroup& group = svData.getDataGroup(0);
    group.add(bamHeader, readsToAdd[0], false, true, true);
    group.add(bamHeader, readsToAdd[3], false, true, true);
    SVSequenceFragmentAssociation association(0, SVEvidenceType::PAIR);
    group.begin().operator*().svLink.push_back(association);

    SVCandidateSetSequenceFragmentSampleGroup& group1 = svData.getDataGroup(1);
    group1.add(bamHeader, readsToAdd[1], false, true, true);
    group1.add(bamHeader, readsToAdd[4], false, true, true);
    group1.begin().operator*().svLink.push_back(association);
    SupportSamples svSupports;
    svSupports.supportSamples.resize(2);
    // Generate SV locus graph file
    // Two nodes are created.
    SVLocus locus1;
    locusAddPair(locus1, 0, 40, 50, 1, 65, 75);
    SVLocus locus2;
    locusAddPair(locus2, 1, 65, 75, 0, 40, 50);
    SVLocusSetOptions sopt;
    sopt.minMergeEdgeObservations = 10;
    SVLocusSet set1(sopt, bamHeader, {bamFileName, bamFileName});
    set1.merge(locus1);
    set1.merge(locus2);
    set1.checkState(true,true);
    TestFilenameMaker testFilenameMaker;
    std::string graphFilename = testFilenameMaker.getFilename();
    const char* testFilenamePtr(graphFilename.c_str());
    // serialize
    set1.save(testFilenamePtr);
    SVCandidateProcessor candidateProcessor(options, scanner.operator*(), programName.c_str(),
                                            version.c_str(), set1, edgeTracker, edgeStatMan);
    candidateProcessor.evaluateCandidates(edgeInfo, mjSvs, svData, svSupports);
    TestSVCandidateProcessor testSVCandidateProcessor;
    testSVCandidateProcessor.flushStreams(candidateProcessor);

    std::ifstream diploidFile(options.diploidOutputFilename);
    std::string line;
    int count(0);
    // Following records are expected:
    // Record-1:    chrFoo	41	MantaBND:0:0:1:0:0:0:0	C	C[chrBar:66[	.	.
    //              SVTYPE=BND;MATEID=MantaBND:0:0:1:0:0:0:1;CIPOS=0,9;HOMLEN=9;HOMSEQ=TCCATGCAT;BND_PAIR_COUNT=0;PAIR_COUNT=1
    // Record-2:    chrBar	66	MantaBND:0:0:1:0:0:0:1	T	]chrFoo:41]T	.	.
    //              SVTYPE=BND;MATEID=MantaBND:0:0:1:0:0:0:0;CIPOS=0,9;HOMLEN=9;HOMSEQ=ATCCTTCTT;BND_PAIR_COUNT=0;PAIR_COUNT=1
    // Record-3:    chrFoo	11	MantaBND:0:0:1:0:0:0:0	C	C]chrBar:110]	0	MinQUAL;NoPairSupport;SampleFT
    //              SVTYPE=BND;MATEID=MantaBND:0:0:1:0:0:0:1;CIPOS=0,9;HOMLEN=9;HOMSEQ=TATCACCCT;BND_DEPTH=3;MATE_BND_DEPTH=3
    //              GT:FT:GQ:PL:PR:SR	0/0:HomRef:48:0,0,0:0,0:0,0
    // Record-4:    chrBar	101	MantaBND:0:0:1:0:0:0:1	G	G]chrFoo:20]	0	MinQUAL;NoPairSupport;SampleFT
    //              SVTYPE=BND;MATEID=MantaBND:0:0:1:0:0:0:0;CIPOS=0,9;HOMLEN=9;HOMSEQ=CAGATGCCA;BND_DEPTH=3;MATE_BND_DEPTH=3
    //              GT:FT:GQ:PL:PR:SR	0/0:HomRef:48:0,0,0:0,0:0,0
    if (diploidFile.is_open())
    {
        while (std::getline(diploidFile, line))
        {
            // ignoring header lines
            if (line.find("#") == std::string::npos)
            {
                count ++;
            }
            if (count >= 1) {
                switch (count) {
                    case 1:
                        BOOST_REQUIRE(line.find("C[chrBar:66[") != std::string::npos);
                        break;
                    case 2:
                        BOOST_REQUIRE(line.find("]chrFoo:41]T") != std::string::npos);
                        break;
                    case 3:
                        BOOST_REQUIRE(line.find("C]chrBar:110]") != std::string::npos);
                        break;
                    case 4:
                        BOOST_REQUIRE(line.find("G]chrFoo:20]") != std::string::npos);
                        break;
                    default:
                        BOOST_THROW_EXCEPTION(illumina::common::GeneralException("less than 1 or more than 4 "
                                                                                 "records are not possible in the VCF"));
                }
            }

        }
        diploidFile.close();
    }
    BOOST_REQUIRE_EQUAL(count, 4);

    std::ifstream somaticFile(options.somaticOutputFilename);
    count = 0;
    // Following records are expected:
    // Record-1:    chrFoo	41	MantaBND:0:0:1:0:0:0:0	C	C[chrBar:66[	.	.
    //              SVTYPE=BND;MATEID=MantaBND:0:0:1:0:0:0:1;CIPOS=0,9;HOMLEN=9;HOMSEQ=TCCATGCAT;BND_PAIR_COUNT=0;PAIR_COUNT=1
    // Record-2:    chrBar	66	MantaBND:0:0:1:0:0:0:1	T	]chrFoo:41]T	.	.
    //              SVTYPE=BND;MATEID=MantaBND:0:0:1:0:0:0:0;CIPOS=0,9;HOMLEN=9;HOMSEQ=ATCCTTCTT;BND_PAIR_COUNT=0;PAIR_COUNT=1
    // Record-3:    chrFoo	11	MantaBND:0:0:1:0:0:0:0	C	C]chrBar:110]	0	MinQUAL;NoPairSupport;SampleFT
    //              SVTYPE=BND;MATEID=MantaBND:0:0:1:0:0:0:1;CIPOS=0,9;HOMLEN=9;HOMSEQ=TATCACCCT;BND_DEPTH=3;MATE_BND_DEPTH=3
    //              GT:FT:GQ:PL:PR:SR	0/0:HomRef:48:0,0,0:0,0:0,0
    // Record-4:    chrBar	101	MantaBND:0:0:1:0:0:0:1	G	G]chrFoo:20]	0	MinQUAL;NoPairSupport;SampleFT
    //              SVTYPE=BND;MATEID=MantaBND:0:0:1:0:0:0:0;CIPOS=0,9;HOMLEN=9;HOMSEQ=CAGATGCCA;BND_DEPTH=3;MATE_BND_DEPTH=3
    //              GT:FT:GQ:PL:PR:SR	0/0:HomRef:48:0,0,0:0,0:0,0
    if (somaticFile.is_open())
    {
        while (std::getline(somaticFile, line))
        {
            if (line.find("#") == std::string::npos)
            {
                count ++;
            }
            if (count >= 1) {
                switch (count) {
                    case 1:
                        BOOST_REQUIRE(line.find("C[chrBar:66[") != std::string::npos);
                        break;
                    case 2:
                        BOOST_REQUIRE(line.find("]chrFoo:41]T") != std::string::npos);
                        break;
                    case 3:
                        BOOST_REQUIRE(line.find("C]chrBar:110]") != std::string::npos);
                        break;
                    case 4:
                        BOOST_REQUIRE(line.find("G]chrFoo:20]") != std::string::npos);
                        break;
                    default:
                        BOOST_THROW_EXCEPTION(illumina::common::GeneralException("less than 1 or more than 4 "
                                                                                 "records are not possible in the VCF"));
                }
            }
        }
        somaticFile.close();
    }
    BOOST_REQUIRE_EQUAL(count, 4);
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()