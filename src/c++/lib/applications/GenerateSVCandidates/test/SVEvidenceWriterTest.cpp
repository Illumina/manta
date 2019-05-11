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

#include "boost/make_unique.hpp"
#include "boost/test/unit_test.hpp"
#include "manta/BamStreamerUtils.hpp"
#include "test/testAlignmentDataUtil.hpp"
#include "test/testFileMakers.hpp"
#include "test/testUtil.hpp"

#include "SVEvidenceWriter.hpp"

BOOST_AUTO_TEST_SUITE(SVEvidenceWriter_test_suite)

// Create Temporary bam streams of a bam file which contains
// two bam record.
struct BamStream {
  BamStream()
  {
    const bam_header_info bamHeader(buildTestBamHeader());

    bam_record bamRecord1;
    buildTestBamRecord(bamRecord1, 0, 200, 1, 300, 200, 15, "200M");
    bamRecord1.toggle_is_first();
    bamRecord1.set_qname("bamRecord1");
    bam_record bamRecord2;
    buildTestBamRecord(bamRecord2, 1, 300, 0, 350, 150, 15, "150M");
    bamRecord2.set_qname("bamRecord2");
    bamRecord2.toggle_is_second();
    readsToAdd.push_back(bamRecord1);
    readsToAdd.push_back(bamRecord2);
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

// check the evidence bam whether datas are as expected. Need to verify the following two fields:
// 1. Read ID
// 2. ZM tag : SV information for that read
// The format of ZM tag should be {SVID1|supportType1|supportType2, SVID2|supportType1|supportType2}
// For example: if a spanning candidate supports deletion, ZM tag will be {DEL_1|PR}
static void checkEvidenceBam(
    const GenomeInterval& genomeInterval,
    const std::string&    bamFileName,
    std::string           expectedZMTagValue,
    std::string           expectedReadName)
{
  // Get the test reference file name
  const std::string referenceFilename = getTestReferenceFilename();
  // List of bam files
  std::vector<std::string>                   bamFilenames = {bamFileName};
  std::vector<std::shared_ptr<bam_streamer>> tempBamStreams;
  openBamStreams(referenceFilename, bamFilenames, tempBamStreams);
  bam_streamer& bamStreamer(*tempBamStreams[0]);

  // Read bam records which are intersecting with the genomeInterval
  bamStreamer.resetRegion(
      genomeInterval.tid, genomeInterval.range.begin_pos(), genomeInterval.range.end_pos());

  // Verify the readID and ZM tag.
  while (bamStreamer.next()) {
    const bam_record* origBamRec(bamStreamer.get_record_ptr());
    BOOST_REQUIRE_EQUAL(origBamRec->qname(), expectedReadName);
    BOOST_REQUIRE_EQUAL(origBamRec->get_string_tag("ZM"), expectedZMTagValue);
  }
}

// Test whether addNewSV api records a single read that support one or more SVs
// for evidence-BAM output correctly or not.
BOOST_AUTO_TEST_CASE(test_AddNewSV)
{
  SVEvidenceWriterRead supportRead;
  supportRead.tid = 0;
  supportRead.pos = 10345;
  // Adding new SV with ID INS_1 and support type PR(spanning pair)
  supportRead.addNewSV("INS_1", "PR");
  // Adding new SV with ID INS_1 and support type SR(split pair)
  supportRead.addNewSV("INS_1", "SR");
  // Adding new SV with ID INS_2 and support type PR(spanning pair)
  supportRead.addNewSV("INS_2", "PR");

  std::map<std::string, std::set<std::string>> svSupportType = supportRead.SVs;
  BOOST_REQUIRE_EQUAL(svSupportType.size(), 2);
  BOOST_REQUIRE_EQUAL(svSupportType["INS_1"].size(), 2);
  BOOST_REQUIRE_EQUAL(svSupportType["INS_2"].size(), 1);
}

// Compare two support reads. Comparison should be based on chromosome number first
// then position.
BOOST_AUTO_TEST_CASE(test_CompareSupportReads)
{
  SVEvidenceWriterRead supportRead1;
  supportRead1.tid = 0;
  supportRead1.pos = 10345;

  SVEvidenceWriterRead supportRead2;
  supportRead2.tid = 1;
  supportRead2.pos = 10345;

  SVEvidenceWriterRead supportRead3;
  supportRead3.tid = 0;
  supportRead3.pos = 10340;
  // Chromosome number of supportRead1 is less than Chromosome number of supportRead2
  BOOST_REQUIRE(supportRead1 < supportRead2);
  BOOST_REQUIRE(!(supportRead2 < supportRead1));
  // Postion of supportRead3 is less than position of supportRead1
  BOOST_REQUIRE(!(supportRead1 < supportRead3));
  BOOST_REQUIRE(supportRead3 < supportRead1);
}

/// Records a single fragment (read1 & read2)
/// that support one or more SVs for evidence-BAM output,
/// indicating the evidence type.
BOOST_AUTO_TEST_CASE(test_SupportFragment)
{
  SVEvidenceWriterReadPair suppFragment1;

  bam_record bamRecord1;
  buildTestBamRecord(bamRecord1, 0, 200, 1, 300, 200, 15, "200M");
  bamRecord1.toggle_is_first();

  // As bamRecord1 is mate1 read, so read1 of supportFragment
  // contains this read information and read2 of supportFragment
  // contains mate information
  suppFragment1.setReads(bamRecord1);
  BOOST_REQUIRE_EQUAL(suppFragment1.read1.tid, 0);
  BOOST_REQUIRE_EQUAL(suppFragment1.read1.pos, 201);
  BOOST_REQUIRE_EQUAL(suppFragment1.read2.tid, 1);
  BOOST_REQUIRE_EQUAL(suppFragment1.read2.pos, 301);

  // Adding spanning read support. It will add "PR" to both
  // read1 and read2.
  suppFragment1.addSpanningSupport("INS_1");
  // Check the size and value after adding spanning support
  BOOST_REQUIRE_EQUAL(suppFragment1.read1.SVs.size(), 1);
  BOOST_REQUIRE_EQUAL(suppFragment1.read2.SVs.size(), 1);
  BOOST_REQUIRE_EQUAL(suppFragment1.read1.SVs["INS_1"].size(), 1);
  BOOST_REQUIRE_EQUAL(suppFragment1.read2.SVs["INS_1"].size(), 1);
  BOOST_REQUIRE_EQUAL(*(suppFragment1.read1.SVs["INS_1"].begin()), "PR");
  BOOST_REQUIRE_EQUAL(*(suppFragment1.read2.SVs["INS_1"].begin()), "PR");

  SVEvidenceWriterReadPair suppFragment2;
  // Adding Split candidate support. If it is mate-1 read, it will add
  // SR to read1 and SRM to read2, if it is mate-2 read, it will add
  // SRM to read1 and SR to read2
  suppFragment2.addSplitSupport(true, "INS_1");  // mate-1 read
  // Check the size after adding split support support to read1
  BOOST_REQUIRE_EQUAL(suppFragment2.read1.SVs.size(), 1);
  BOOST_REQUIRE_EQUAL(suppFragment2.read2.SVs.size(), 1);
  BOOST_REQUIRE_EQUAL(suppFragment2.read1.SVs["INS_1"].size(), 1);
  BOOST_REQUIRE_EQUAL(suppFragment2.read2.SVs["INS_1"].size(), 1);
  BOOST_REQUIRE_EQUAL(*(suppFragment2.read1.SVs["INS_1"].begin()), "SR");
  BOOST_REQUIRE_EQUAL(*(suppFragment2.read2.SVs["INS_1"].begin()), "SRM");

  SVEvidenceWriterReadPair suppFragment3;
  suppFragment3.addSplitSupport(false, "INS_1");  // mate-2 read
  BOOST_REQUIRE_EQUAL(suppFragment3.read1.SVs.size(), 1);
  BOOST_REQUIRE_EQUAL(suppFragment3.read2.SVs.size(), 1);
  BOOST_REQUIRE_EQUAL(suppFragment3.read1.SVs["INS_1"].size(), 1);
  BOOST_REQUIRE_EQUAL(suppFragment3.read2.SVs["INS_1"].size(), 1);
  BOOST_REQUIRE_EQUAL(*(suppFragment3.read1.SVs["INS_1"].begin()), "SRM");
  BOOST_REQUIRE_EQUAL(*(suppFragment3.read2.SVs["INS_1"].begin()), "SR");
}

// Test SVEvidenceWriterSampleData which records all supporting fragments
// that support one or more SVs for evidence-BAM output.
BOOST_AUTO_TEST_CASE(test_SupportFragments)
{
  SVEvidenceWriterSampleData suppFragments;
  bam_record                 bamRecord1;
  buildTestBamRecord(bamRecord1, 0, 200, 1, 300, 200, 15, "200M");
  bamRecord1.toggle_is_first();
  bamRecord1.set_qname("bamRecord1");

  // As bamRecord1 is mate-1 read, read1 should contain this read's information
  // and read2 should contain it's mate information.
  SVEvidenceWriterReadPair suppFragment1 = suppFragments.getSupportFragment(bamRecord1);
  BOOST_REQUIRE_EQUAL(suppFragment1.read1.tid, 0);
  BOOST_REQUIRE_EQUAL(suppFragment1.read1.pos, 201);
  BOOST_REQUIRE_EQUAL(suppFragment1.read2.tid, 1);
  BOOST_REQUIRE_EQUAL(suppFragment1.read2.pos, 301);

  bam_record bamRecord2;
  buildTestBamRecord(bamRecord2, 1, 350, 0, 250, 100, 15, "100M");
  bamRecord2.toggle_is_second();
  bamRecord2.set_qname("bamRecord2");

  SVCandidateSetSequenceFragment candidateSetSequenceFragment1;
  candidateSetSequenceFragment1.read1.bamrec = bamRecord1;

  // Testing same information as above using SVCandidateSetSequenceFragment
  SVEvidenceWriterReadPair suppFragment2 = suppFragments.getSupportFragment(candidateSetSequenceFragment1);
  BOOST_REQUIRE_EQUAL(suppFragment2.read1.tid, 0);
  BOOST_REQUIRE_EQUAL(suppFragment2.read1.pos, 201);
  BOOST_REQUIRE_EQUAL(suppFragment2.read2.tid, 1);
  BOOST_REQUIRE_EQUAL(suppFragment2.read2.pos, 301);

  SVCandidateSetSequenceFragment candidateSetSequenceFragment2;
  candidateSetSequenceFragment2.read2.bamrec = bamRecord2;
  SVEvidenceWriterReadPair suppFragment3 = suppFragments.getSupportFragment(candidateSetSequenceFragment2);
  // As bamRecord2 is mate-2 read, read1 should contain it's mate information
  // and read2 should contain this read's information.
  BOOST_REQUIRE_EQUAL(suppFragment3.read1.tid, 0);
  BOOST_REQUIRE_EQUAL(suppFragment3.read1.pos, 251);
  BOOST_REQUIRE_EQUAL(suppFragment3.read2.tid, 1);
  BOOST_REQUIRE_EQUAL(suppFragment3.read2.pos, 351);
}

// Test SVEvidenceWriterData which is vector of support fragments
BOOST_AUTO_TEST_CASE(test_SupportSamples)
{
  SVEvidenceWriterData       svEvidenceWriterData(2);
  SVEvidenceWriterSampleData suppFragments1;
  bam_record                 bamRecord1;
  buildTestBamRecord(bamRecord1, 0, 200, 1, 300, 200, 15, "200M");
  bamRecord1.toggle_is_first();
  bamRecord1.set_qname("bamRecord1");
  svEvidenceWriterData.sampleData[0] = suppFragments1;

  SVEvidenceWriterSampleData suppFragments2;
  bam_record                 bamRecord2;
  buildTestBamRecord(bamRecord2, 1, 1200, 0, 1300, 200, 15, "200M");
  bamRecord2.toggle_is_first();
  bamRecord2.set_qname("bamRecord2");
  svEvidenceWriterData.sampleData[1] = suppFragments2;

  // Check information for bamRecord1
  SVEvidenceWriterReadPair suppFragment1 =
      svEvidenceWriterData.getSampleData(0).getSupportFragment(bamRecord1);
  BOOST_REQUIRE_EQUAL(suppFragment1.read1.tid, 0);
  BOOST_REQUIRE_EQUAL(suppFragment1.read1.pos, 201);
  BOOST_REQUIRE_EQUAL(suppFragment1.read2.tid, 1);
  BOOST_REQUIRE_EQUAL(suppFragment1.read2.pos, 301);

  // Check information for bamRecord2
  SVEvidenceWriterReadPair suppFragment2 =
      svEvidenceWriterData.getSampleData(1).getSupportFragment(bamRecord2);
  BOOST_REQUIRE_EQUAL(suppFragment2.read1.tid, 1);
  BOOST_REQUIRE_EQUAL(suppFragment2.read1.pos, 1201);
  BOOST_REQUIRE_EQUAL(suppFragment2.read2.tid, 0);
  BOOST_REQUIRE_EQUAL(suppFragment2.read2.pos, 1301);
}

// Take the evidence read from the original bam and add new ZM (SV information)tag to the bam record.
// check the evidence bam whether datas are as expected. Need to verify the following two fields:
// 1. Read ID
// 2. ZM tag : SV information for that read
// The format of ZM tag should be {SVID1|supportType1|supportType2, SVID2|supportType1|supportType2}
// For example: if a spanning candidate supports deletion, ZM tag will be {DEL_1|PR}
BOOST_AUTO_TEST_CASE(test_ProcessRecords)
{
  GenomeInterval               genomeInterval1(0, 100, 220);
  GenomeInterval               genomeInterval2(1, 300, 350);
  const bam_header_info        bamHeader(buildTestBamHeader());
  const HtslibBamHeaderManager bamHeaderManager(bamHeader.chrom_data);
  const BamFilenameMaker       bamFileNameMaker;
  const std::string&           bamFileName(bamFileNameMaker.getFilename());

  SVEvidenceWriterSampleData suppFragments;
  SVEvidenceWriterReadPair   suppFragment1;
  suppFragment1.setReads(readsToAdd[0]);
  suppFragment1.addSpanningSupport("INS_1");
  SVEvidenceWriterReadPair suppFragment2;
  suppFragment2.setReads(readsToAdd[1]);
  suppFragment2.addSpanningSupport("DEL_1");
  suppFragments.supportFrags[readsToAdd[0].qname()] = suppFragment1;
  suppFragments.supportFrags[readsToAdd[1].qname()] = suppFragment2;

  {
    // this code block forces bamWriter to flush at block end:
    SynchronizedBamWriter bamWriter(bamFileName.c_str(), bamHeaderManager.get());
    // Process the bam record and add ZM tag for all reads which are intersection to genomeInterval1
    SVEvidenceWriter::processBamRecords(
        bamStream.operator*(), genomeInterval1, suppFragments.supportFrags, bamWriter);
  }

  // create bam index
  const int indexStatus1 = bam_index_build(bamFileName.c_str(), 0);
  // check whether .bai file for bam file is build successfully or not.
  BOOST_REQUIRE(indexStatus1 >= 0);

  // check the evidence bam for bamRecord1 as genomeInterval1 intersects with bamRecord1.
  checkEvidenceBam(genomeInterval1, bamFileName, "INS_1|PR", "bamRecord1");

  {
    // this code block forces bamWriter to flush at block end:
    SynchronizedBamWriter bamWriter(bamFileName.c_str(), bamHeaderManager.get());
    // Process the bam record and add ZM tag for all reads which are intersection to genomeInterval2
    SVEvidenceWriter::processBamRecords(
        bamStream.operator*(), genomeInterval1, suppFragments.supportFrags, bamWriter);
  }

  // create bam index
  const int indexStatus2 = bam_index_build(bamFileName.c_str(), 0);
  // check whether .bai file for bam file is build successfully or not.
  BOOST_REQUIRE(indexStatus2 >= 0);

  // check the evidence bam for bamRecord2 as genomeInterval2 intersects with bamRecord2.
  checkEvidenceBam(genomeInterval2, bamFileName, "DEL_1|PR", "bamRecord2");
}

// Test the Writing of evidence bam.
BOOST_AUTO_TEST_CASE(test_writeEvidenceBam)
{
  GenomeInterval genomeInterval1(0, 100, 220);
  GenomeInterval genomeInterval2(1, 300, 350);

  const bam_header_info        bamHeader(buildTestBamHeader());
  const HtslibBamHeaderManager bamHeaderManager(bamHeader.chrom_data);
  const BamFilenameMaker       bamFileNameMaker;
  const std::string&           bamFileName(bamFileNameMaker.getFilename());

  SVEvidenceWriterSampleData suppFragments;
  SVEvidenceWriterReadPair   suppFragment1;
  suppFragment1.setReads(readsToAdd[0]);
  suppFragment1.addSpanningSupport("INS_1");
  SVEvidenceWriterReadPair suppFragment2;
  suppFragment2.setReads(readsToAdd[1]);
  suppFragment2.addSpanningSupport("DEL_1");

  suppFragments.supportFrags[readsToAdd[0].qname()] = suppFragment1;
  suppFragments.supportFrags[readsToAdd[1].qname()] = suppFragment2;

  {
    // this code block forces bamWriter to flush at block end:
    SynchronizedBamWriter bamWriter(bamFileName.c_str(), bamHeaderManager.get());
    SVEvidenceWriter::writeSupportBam(bamStream, suppFragments, bamWriter);
  }

  // create bam index
  const int indexStatus = bam_index_build(bamFileName.c_str(), 0);
  // check whether .bai file for bam file is build successfully or not.
  BOOST_REQUIRE(indexStatus >= 0);

  // check the evidence bam as mentioned in the doc in test_ProcessRecords
  checkEvidenceBam(genomeInterval1, bamFileName, "INS_1|PR", "bamRecord1");
  checkEvidenceBam(genomeInterval2, bamFileName, "DEL_1|PR", "bamRecord2");
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
