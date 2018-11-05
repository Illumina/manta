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

#include "boost/test/unit_test.hpp"
#include "boost/make_unique.hpp"
#include "test/testAlignmentDataUtil.hh"
#include "test/testAlignmentDataUtil.hh"
#include "test/testFileMakers.hh"
#include "test/testUtil.hh"
#include "manta/BamStreamerUtils.hh"

#include "SVSupports.hh"

BOOST_AUTO_TEST_SUITE( SVSupports_test_suite )

// Create Temporary bam streams of a bam file which contains
// two bam record.
struct OpenBamStream
{
    OpenBamStream()
    {
        const bam_header_info bamHeader(buildTestBamHeader());

        bam_record bam_record1;
        buildTestBamRecord(bam_record1, 0, 200, 1, 300, 200, 15, "200M");
        bam_record1.toggle_is_first();
        bam_record1.set_qname("bam_record1");
        bam_record bam_record2;
        buildTestBamRecord(bam_record2, 1, 300, 0, 350, 150, 15, "150M");
        bam_record2.set_qname("bam_record2");
        bam_record2.toggle_is_second();
        readsToAdd.push_back(bam_record1);
        readsToAdd.push_back(bam_record2);
        buildTestBamFile(bamHeader, readsToAdd, _bamFilename());

        const std::string referenceFilename = getTestReferenceFilename();
        std::vector<std::string> bamFilenames = { _bamFilename() };
        openBamStreams(referenceFilename, bamFilenames, bamStreams);
    }

    std::vector<std::shared_ptr<bam_streamer>> bamStreams;
    std::vector<bam_record> readsToAdd;
    support_fragments_t fragmentsMap;
private:

    const std::string&
    _bamFilename() const
    {
        return _bamFilenameMaker.getFilename();
    }
    const BamFilenameMaker _bamFilenameMaker;
};

BOOST_FIXTURE_TEST_SUITE( SVSupports_test_suite, OpenBamStream)


// check the evidence bam whether datas are as expected. Need to verify the following two fields:
// 1. Read ID
// 2. ZM tag : SV information for that read
// The format of ZM tag should be {SVID1|supportType1|supportType2, SVID2|supportType1|supportType2}
// For example: if a spanning candidate supports deletion, ZM tag will be {DEL|PR}
static void checkEvidenceBam(const GenomeInterval &genomeInterval,
                             const std::string &bamFileName,
                             std::string expectedZMTagValue, std::string expectedReadName)
{
    // Get the test reference file name
    const std::string referenceFilename = getTestReferenceFilename();
    // List of bam files
    std::vector<std::string> bamFilenames = { bamFileName};
    std::vector<std::shared_ptr<bam_streamer>> tempBamStreams;
    openBamStreams(referenceFilename, bamFilenames, tempBamStreams);
    bam_streamer& bamStreamer(*tempBamStreams[0]);

    // Read bam records which are intersecting with the genomeInterval
    bamStreamer.resetRegion(genomeInterval.tid, genomeInterval.range.begin_pos(), genomeInterval.range.end_pos());

    // Verify the readID and ZM tag.
    while (bamStreamer.next())
    {
        const bam_record* origBamRec(bamStreamer.get_record_ptr());
        bam_record bamRec(*origBamRec);
        BOOST_REQUIRE_EQUAL(bamRec.qname(), expectedReadName);
        BOOST_REQUIRE_EQUAL(bamRec.get_string_tag("ZM"), expectedZMTagValue);
    }
}

// Test whether addNewSV api records a single read that support one or more SVs
// for evidence-BAM output correctly or not.
BOOST_AUTO_TEST_CASE( test_AddNewSV)
{
    SupportRead supportRead;
    supportRead.tid = 0;
    supportRead.pos = 10345;
    // Adding new SV with ID INS and support type PR(spanning pair)
    supportRead.addNewSV("INS", "PR");
    // Adding new SV with ID INS and support type SR(split pair)
    supportRead.addNewSV("INS", "SR");
    // Adding new SV with ID DEL and support type PR(spanning pair)
    supportRead.addNewSV("DEL", "PR");

    std::map<std::string, std::set<std::string>> sv_supportType = supportRead.SVs;
    BOOST_REQUIRE_EQUAL(sv_supportType.size(), 2);
    BOOST_REQUIRE_EQUAL(sv_supportType["INS"].size(), 2);
    BOOST_REQUIRE_EQUAL(sv_supportType["DEL"].size(), 1);
}

// Compare two support reads. Comparison should be based on chromosome number first
// then position.
BOOST_AUTO_TEST_CASE( test_CompareSupportReads)
{
    SupportRead supportRead1;
    supportRead1.tid = 0;
    supportRead1.pos = 10345;

    SupportRead supportRead2;
    supportRead2.tid = 1;
    supportRead2.pos = 10345;

    SupportRead supportRead3;
    supportRead3.tid = 0;
    supportRead3.pos = 10340;
    // Chromosome number of supportRead1 is less than Chromosome number of supportRead2
    BOOST_REQUIRE(supportRead1 < supportRead2);
    BOOST_REQUIRE(!(supportRead2 < supportRead1));
    // Postion of supportRead3 is less than position of supportRead1
    BOOST_REQUIRE(!(supportRead1 < supportRead3));
    BOOST_REQUIRE(supportRead3 < supportRead1);
}

// test the output stream
BOOST_AUTO_TEST_CASE( test_PrintSupportRead)
{
    SupportRead supportRead;
    supportRead.tid = 0;
    supportRead.pos = 10345;
    supportRead.addNewSV("INS", "PR");
    supportRead.addNewSV("INS", "SR");
    supportRead.addNewSV("DEL", "PR");
    std::cout << supportRead;
}

/// Records a single fragment (read1 & read2)
/// that support one or more SVs for evidence-BAM output,
/// indicating the evidence type.
BOOST_AUTO_TEST_CASE( test_SupportFragment )
{
    SupportFragment supportFragment;

    bam_record bam_record1;
    buildTestBamRecord(bam_record1, 0, 200, 1, 300, 200, 15, "200M");
    bam_record1.toggle_is_first();

    // As bam_record1 is mate1 read, so read1 of supportFragment
    // contains this read information and read2 of supportFragment
    // contains mate information
    supportFragment.setReads(bam_record1);
    BOOST_REQUIRE_EQUAL(supportFragment.read1.tid, 0);
    BOOST_REQUIRE_EQUAL(supportFragment.read1.pos, 201);
    BOOST_REQUIRE_EQUAL(supportFragment.read2.tid, 1);
    BOOST_REQUIRE_EQUAL(supportFragment.read2.pos, 301);

    // Adding spanning read support. It will add "PR" to both
    // read1 and read2.
    supportFragment.addSpanningSupport("INS");
    supportFragment.addSpanningSupport("DEL");

    // Check the size after adding spanning support
    BOOST_REQUIRE_EQUAL(supportFragment.read1.SVs.size(), 2);
    BOOST_REQUIRE_EQUAL(supportFragment.read2.SVs.size(), 2);
    BOOST_REQUIRE_EQUAL(supportFragment.read1.SVs["INS"].size(), 1);
    BOOST_REQUIRE_EQUAL(supportFragment.read1.SVs["DEL"].size(), 1);
    BOOST_REQUIRE_EQUAL(supportFragment.read2.SVs["INS"].size(), 1);
    BOOST_REQUIRE_EQUAL(supportFragment.read2.SVs["DEL"].size(), 1);

    // Adding Split candidate support. If it is mate-1 read, it will add
    // SR to read1 and SRM to read2, if it is mate-2 read, it will add
    // SRM to read1 and SR to read2
    supportFragment.addSplitSupport(true, "INS");
    // Check the size after adding split support support to read1
    BOOST_REQUIRE_EQUAL(supportFragment.read1.SVs.size(), 2);
    BOOST_REQUIRE_EQUAL(supportFragment.read2.SVs.size(), 2);
    BOOST_REQUIRE_EQUAL(supportFragment.read1.SVs["INS"].size(), 2);
    BOOST_REQUIRE_EQUAL(supportFragment.read1.SVs["DEL"].size(), 1);
    BOOST_REQUIRE_EQUAL(supportFragment.read2.SVs["INS"].size(), 2);
    BOOST_REQUIRE_EQUAL(supportFragment.read2.SVs["DEL"].size(), 1);

    supportFragment.addSplitSupport(false, "INS");
    // Check the size after adding split support support to read2
    BOOST_REQUIRE_EQUAL(supportFragment.read1.SVs.size(), 2);
    BOOST_REQUIRE_EQUAL(supportFragment.read2.SVs.size(), 2);
    BOOST_REQUIRE_EQUAL(supportFragment.read1.SVs["INS"].size(), 3);
    BOOST_REQUIRE_EQUAL(supportFragment.read1.SVs["DEL"].size(), 1);
    BOOST_REQUIRE_EQUAL(supportFragment.read2.SVs["INS"].size(), 3);
    BOOST_REQUIRE_EQUAL(supportFragment.read2.SVs["DEL"].size(), 1);

    // check the output stream
    std::cout << supportFragment << std::endl;
}

// Test SupportFragments which records all supporting fragments
// that support one or more SVs for evidence-BAM output.
BOOST_AUTO_TEST_CASE( test_SupportFragments )
{
    SupportFragments supportFragments;
    bam_record bam_record1;
    buildTestBamRecord(bam_record1, 0, 200, 1, 300, 200, 15, "200M");
    bam_record1.toggle_is_first();
    bam_record1.set_qname("bam_record_1");

    // As bam_record1 is mate-1 read, read1 should contain this read's information
    // and read2 should contain it's mate information.
    SupportFragment supportFragment1 = supportFragments.getSupportFragment(bam_record1);
    BOOST_REQUIRE_EQUAL(supportFragment1.read1.tid, 0);
    BOOST_REQUIRE_EQUAL(supportFragment1.read1.pos, 201);
    BOOST_REQUIRE_EQUAL(supportFragment1.read2.tid, 1);
    BOOST_REQUIRE_EQUAL(supportFragment1.read2.pos, 301);

    bam_record bam_record2;
    buildTestBamRecord(bam_record2, 1, 350, 0, 250, 200, 15, "100M");
    bam_record2.toggle_is_second();
    bam_record2.set_qname("bam_record_2");

    SVCandidateSetSequenceFragment svCandidateSetSequenceFragment1;
    svCandidateSetSequenceFragment1.read1.bamrec = bam_record1;

    // Testing same information as above using SVCandidateSetSequenceFragment
    SupportFragment supportFragment2 = supportFragments.getSupportFragment(svCandidateSetSequenceFragment1);
    BOOST_REQUIRE_EQUAL(supportFragment2.read1.tid, 0);
    BOOST_REQUIRE_EQUAL(supportFragment2.read1.pos, 201);
    BOOST_REQUIRE_EQUAL(supportFragment2.read2.tid, 1);
    BOOST_REQUIRE_EQUAL(supportFragment2.read2.pos, 301);

    SVCandidateSetSequenceFragment svCandidateSetSequenceFragment2;
    svCandidateSetSequenceFragment2.read2.bamrec = bam_record2;
    SupportFragment supportFragment3 = supportFragments.getSupportFragment(svCandidateSetSequenceFragment2);
    // As bam_record2 is mate-2 read, read1 should contain it's mate information
    // and read2 should contain this read's information.
    BOOST_REQUIRE_EQUAL(supportFragment3.read1.tid, 0);
    BOOST_REQUIRE_EQUAL(supportFragment3.read1.pos, 251);
    BOOST_REQUIRE_EQUAL(supportFragment3.read2.tid, 1);
    BOOST_REQUIRE_EQUAL(supportFragment3.read2.pos, 351);

    // check the output stream
    std::cout << supportFragments << std::endl;

}

// Test SupportSamples which is vector of support fragments
BOOST_AUTO_TEST_CASE( test_SupportSamples )
{
    SupportSamples supportSamples;
    supportSamples.supportSamples.resize(1);
    SupportFragments supportFragments1;
    bam_record bam_record1;
    buildTestBamRecord(bam_record1, 0, 200, 1, 300, 200, 15, "200M");
    bam_record1.toggle_is_first();
    bam_record1.set_qname("bam_record1");
    supportSamples.supportSamples[0] = supportFragments1;

    SupportFragments supportFragments2;
    bam_record bam_record2;
    buildTestBamRecord(bam_record2, 1, 1200, 0, 1300, 200, 15, "200M");
    bam_record2.toggle_is_first();
    bam_record2.set_qname("bam_record2");
    supportSamples.supportSamples[0] = supportFragments2;

    // Check information for bam_record1
    SupportFragment supportFragment1 = supportSamples.getSupportFragments(0).getSupportFragment(bam_record1);
    BOOST_REQUIRE_EQUAL(supportFragment1.read1.tid, 0);
    BOOST_REQUIRE_EQUAL(supportFragment1.read1.pos, 201);
    BOOST_REQUIRE_EQUAL(supportFragment1.read2.tid, 1);
    BOOST_REQUIRE_EQUAL(supportFragment1.read2.pos, 301);

    // Check information for bam_record2
    SupportFragment supportFragment2 = supportSamples.getSupportFragments(0).getSupportFragment(bam_record2);
    BOOST_REQUIRE_EQUAL(supportFragment2.read1.tid, 1);
    BOOST_REQUIRE_EQUAL(supportFragment2.read1.pos, 1201);
    BOOST_REQUIRE_EQUAL(supportFragment2.read2.tid, 0);
    BOOST_REQUIRE_EQUAL(supportFragment2.read2.pos, 1301);

    // check the output stream
    std::cout << supportSamples << std::endl;
}

// Take the evidence read from the original bam and add new ZM (SV information)tag to the bam record.
// check the evidence bam whether datas are as expected. Need to verify the following two fields:
// 1. Read ID
// 2. ZM tag : SV information for that read
// The format of ZM tag should be {SVID1|supportType1|supportType2, SVID2|supportType1|supportType2}
// For example: if a spanning candidate supports deletion, ZM tag will be {DEL|PR}
BOOST_AUTO_TEST_CASE( test_ProcessRecords )
{

    GenomeInterval genomeInterval1(0, 100, 220);
    GenomeInterval genomeInterval2(1, 300, 350);
    const bam_header_info bamHeader(buildTestBamHeader());
    const HtslibBamHeaderManager bamHeaderManager(bamHeader.chrom_data);
    const std::shared_ptr<BamFilenameMaker> bamFileNameMaker(boost::make_unique<BamFilenameMaker>());
    std::string bamFileName = bamFileNameMaker.get()->getFilename();

    SupportFragments fragments;
    SupportFragment supportFragment1;
    supportFragment1.setReads(readsToAdd[0]);
    supportFragment1.addSpanningSupport("INS");
    supportFragment1.addSpanningSupport("DEL");
    SupportFragment supportFragment2;
    supportFragment2.setReads(readsToAdd[1]);
    supportFragment2.addSpanningSupport("DEL");
    fragments.supportFrags[readsToAdd[0].qname()] = supportFragment1;
    fragments.supportFrags[readsToAdd[1].qname()] = supportFragment2;

    bam_dumper bamDumper1(bamFileName.c_str(), bamHeaderManager.get());
    // Process the bam record and add ZM tag for all reads which are intersection to genomeInterval1
    processBamRecords(bamStreams[0].operator*(), genomeInterval1, fragments.supportFrags, bamDumper1);
    bamDumper1.close();
    // creating bam index
    const int indexStatus1 = bam_index_build(bamFileName.c_str(), 0);
    std::cout << "Bam Index is" <<  indexStatus1 << std::endl;
    // check the evidence bam for bam_record1 as genomeInterval1 intersects with bam_record1.
    checkEvidenceBam(genomeInterval1, bamFileName, "DEL|PR,INS|PR", "bam_record1");

    bam_dumper bamDumper2(bamFileName.c_str(), bamHeaderManager.get());
    // Process the bam record and add ZM tag for all reads which are intersection to genomeInterval2
    processBamRecords(bamStreams[0].operator*(), genomeInterval1, fragments.supportFrags, bamDumper2);
    bamDumper2.close();
    const int indexStatus2 = bam_index_build(bamFileName.c_str(), 0);
    std::cout << "Bam Index is" <<  indexStatus2 << std::endl;
    // check the evidence bam for bam_record2 as genomeInterval2 intersects with bam_record2.
    checkEvidenceBam(genomeInterval2, bamFileName, "DEL|PR", "bam_record2");

    // cleanup
    remove(bamFileName.c_str());
    remove((bamFileName + ".bai").c_str());
}

// Test the Writing of evidence bam.
BOOST_AUTO_TEST_CASE( test_writeSupportBam )
{
    GenomeInterval genomeInterval1(0, 100, 220);
    GenomeInterval genomeInterval2(1, 300, 350);

    const bam_header_info bamHeader(buildTestBamHeader());
    const HtslibBamHeaderManager bamHeaderManager(bamHeader.chrom_data);
    const std::shared_ptr<BamFilenameMaker> bamFileNameMaker(boost::make_unique<BamFilenameMaker>());
    std::string bamFileName = bamFileNameMaker.get()->getFilename();

    SupportFragments fragments;
    SupportFragment supportFragment1;
    supportFragment1.setReads(readsToAdd[0]);
    supportFragment1.addSpanningSupport("INS");
    supportFragment1.addSpanningSupport("DEL");
    SupportFragment supportFragment2;
    supportFragment2.setReads(readsToAdd[1]);
    supportFragment2.addSpanningSupport("DEL");

    fragments.supportFrags[readsToAdd[0].qname()] = supportFragment1;
    fragments.supportFrags[readsToAdd[1].qname()] = supportFragment2;

    bam_dumper_ptr bamDumperPtr(new bam_dumper(bamFileName.c_str(), bamHeaderManager.get()));
    // Write both the bam_record1 and bam_record2 in the evidence bam
    writeSupportBam(bamStreams[0], fragments, bamDumperPtr);
    bamDumperPtr.get()->close();
    // creating bam index
    const int indexStatus = bam_index_build(bamFileName.c_str(), 0);
    std::cout << "Bam Index is" <<  indexStatus << std::endl;

    // check the evidence bam as mentioned in the doc in test_ProcessRecords
    checkEvidenceBam(genomeInterval1, bamFileName, "DEL|PR,INS|PR", "bam_record1");
    checkEvidenceBam(genomeInterval2, bamFileName, "DEL|PR", "bam_record2");

    // cleanup
    remove(bamFileName.c_str());
    remove((bamFileName + ".bai").c_str());
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()