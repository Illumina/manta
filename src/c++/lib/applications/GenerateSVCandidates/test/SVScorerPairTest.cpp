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

#include "boost/make_unique.hpp"
#include "boost/test/unit_test.hpp"
#include "manta/BamStreamerUtils.hpp"
#include "test/testAlignmentDataUtil.hpp"
#include "test/testFileMakers.hpp"
#include "test/testSVLocusScanner.hpp"

#include "SVScorer.hpp"
#include "SVScorerPair.cpp"

/// TestSVScorer is a friend of SVScorer. So that can access private
/// methods of SVScorer.
struct TestSVScorer {
  void processExistingAltPairInfo(
      SVScorer&             scorer,
      PairOptions&          pairOptions,
      SVCandidateSetData&   candidateSetData,
      SVCandidate&          candidate,
      SVId&                 id,
      SVEvidence&           evidence,
      SVEvidenceWriterData& suppSamples)
  {
    scorer.processExistingAltPairInfo(pairOptions, candidateSetData, candidate, id, evidence, suppSamples);
  }

  void getSVPairSupport(
      SVScorer&               scorer,
      SVCandidateSetData&     candidateSetData,
      SVCandidateAssemblyData assemblyData,
      SVCandidate&            candidate,
      SVId&                   id,
      SVEvidence&             evidence,
      SVEvidenceWriterData&   suppSamples)
  {
    scorer.getSVPairSupport(candidateSetData, assemblyData, candidate, id, evidence, suppSamples);
  }
};

std::unique_ptr<SVLocusScanner> buildSVLocusScanner(bam_header_info bamHeaderInfo);
// Creating locus scanner for two bam files
std::unique_ptr<SVLocusScanner> buildSVLocusScanner(bam_header_info bamHeaderInfo)
{
  TestFilenameMaker filenameMaker;
  std::string       tmpFileName = filenameMaker.getFilename();
  ReadGroupStatsSet rstats;
  for (unsigned j(0); j < 2; j++) {
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
  ReadScannerOptions opts      = ReadScannerOptions();
  opts.minCandidateVariantSize = 8;
  const ReadScannerOptions&      constRefOpts(opts);
  TestAlignHeaderFileMaker       alignFile1(bamHeaderInfo);
  TestAlignHeaderFileMaker       alignFile2(bamHeaderInfo);
  const std::vector<std::string> alignFilenameVector = {alignFile1.getFilename(), alignFile2.getFilename()};
  return boost::make_unique<SVLocusScanner>(constRefOpts, tmpFileName, alignFilenameVector, false);
}

BOOST_AUTO_TEST_SUITE(SVScorerPair_test_suite)

// Create Temporary bam streams of a bam file which contains
// three bam records.
struct BamStream {
  BamStream()
  {
    const bam_header_info bamHeader(buildTestBamHeader());

    std::string querySeq = "TCTATCACCCATTTTACCACTCACGGGAGCTCTCC";
    // Creating 1st bam
    bam_record bamRecord1;
    buildTestBamRecord(bamRecord1, 0, 9, 0, 70, 35, 15, "35M", querySeq, 125);
    bamRecord1.set_qname("bamRecord1");
    bam_record bamRecord2;
    buildTestBamRecord(bamRecord2, 0, 109, 0, 170, 35, 15, "35M", querySeq, 130);
    bamRecord2.set_qname("bamRecord2");
    bam_record bamRecord3;
    buildTestBamRecord(bamRecord3, 0, 109, 0, 170, 35, 15, "35M", querySeq, 130);
    bamRecord3.set_qname("bamRecord3");
    bam_record bamRecord4;
    buildTestBamRecord(bamRecord4, 0, 109, 0, 175, 35, 15, "35M", querySeq, 100);
    bamRecord4.set_qname("bamRecord4");
    readsToAddForFirstBam.push_back(bamRecord1);
    readsToAddForFirstBam.push_back(bamRecord2);
    readsToAddForFirstBam.push_back(bamRecord3);
    readsToAddForFirstBam.push_back(bamRecord4);
    bamFileName1 = _bamFilename1();
    buildTestBamFile(bamHeader, readsToAddForFirstBam, bamFileName1);

    // Creating 2nd bam
    bam_record bamRecord5;
    buildTestBamRecord(bamRecord5, 0, 10, 0, 70, 35, 15, "35M", querySeq, 125);
    bamRecord5.set_qname("bamRecord5");
    bam_record bamRecord6;
    buildTestBamRecord(bamRecord6, 0, 108, 0, 170, 35, 15, "35M", querySeq, 130);
    bamRecord6.set_qname("bamRecord6");
    bam_record bamRecord7;
    buildTestBamRecord(bamRecord7, 0, 108, 0, 170, 35, 15, "35M", querySeq, 130);
    bamRecord7.set_qname("bamRecord7");
    bam_record bamRecord8;
    buildTestBamRecord(bamRecord8, 0, 109, 0, 175, 35, 15, "35M", querySeq, 100);
    bamRecord8.set_qname("bamRecord8");
    readsToAddForSecondBam.push_back(bamRecord5);
    readsToAddForSecondBam.push_back(bamRecord6);
    readsToAddForSecondBam.push_back(bamRecord7);
    readsToAddForSecondBam.push_back(bamRecord8);
    bamFileName2 = _bamFilename2();
    buildTestBamFile(bamHeader, readsToAddForSecondBam, bamFileName2);
    const std::string                          referenceFilename = getTestReferenceFilename();
    std::vector<std::string>                   bamFilenames      = {bamFileName1, bamFileName2};
    std::vector<std::shared_ptr<bam_streamer>> bamStreams;
    openBamStreams(referenceFilename, bamFilenames, bamStreams);
    bamStream1 = bamStreams[0];
    bamStream2 = bamStreams[1];
  }

  std::shared_ptr<bam_streamer> bamStream1;
  std::vector<bam_record>       readsToAddForFirstBam;
  std::string                   bamFileName1;
  std::shared_ptr<bam_streamer> bamStream2;
  std::vector<bam_record>       readsToAddForSecondBam;
  std::string                   bamFileName2;

private:
  const std::string& _bamFilename1() const { return _bamFilenameMaker1.getFilename(); }

  const std::string&     _bamFilename2() const { return _bamFilenameMaker2.getFilename(); }
  const BamFilenameMaker _bamFilenameMaker1;
  const BamFilenameMaker _bamFilenameMaker2;
};

BOOST_FIXTURE_TEST_SUITE(SVScorerPair_test_suite, BamStream)

// getFragInfo api fills read information as accurate as possible based
// on either only mate-1 or mate-2 or both. Following cases need to be tested
// 1. When only mate-1 read available.
// 2. When mate-1 and mate-2 both are available.
BOOST_AUTO_TEST_CASE(test_getFragInfo)
{
  SVCandidateSetSequenceFragment fragment1;
  SpanReadInfo                   read1;
  SpanReadInfo                   read2;

  std::string readSeq1 = "GGTATTTTCGTCTGGGGGGT";
  std::string readSeq2 = "GGTATTTTCG";
  // when only mate-1 read is available.
  bam_record bamRecord1;
  buildTestBamRecord(bamRecord1, 0, 200, 0, 270, 20, 15, "20M", readSeq1);
  bamRecord1.set_qname("Read-1");
  fragment1.read1.bamrec = bamRecord1;
  getFragInfo(fragment1, read1, read2);
  BOOST_REQUIRE(read1.isFwdStrand);
  BOOST_REQUIRE_EQUAL(read1.readSize, 20);
  BOOST_REQUIRE_EQUAL(read1.interval, GenomeInterval(0, 200, 220));
  // approximate values of mate-2 as mate-2 is not available
  BOOST_REQUIRE(!read2.isFwdStrand);
  BOOST_REQUIRE_EQUAL(read2.readSize, 20);
  BOOST_REQUIRE_EQUAL(read2.interval, GenomeInterval(0, 270, 290));

  // When mate-1 and mate-2 both are available
  bam_record bamRecord2;
  buildTestBamRecord(bamRecord2, 0, 270, 0, 200, 10, 15, "10M", readSeq2);
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

  // When only mate-2 read is available
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
BOOST_AUTO_TEST_CASE(test_getTerminal)
{
  SpanReadInfo readInfo;
  SpanTerminal terminal;
  // case of forward strand
  readInfo.isFwdStrand = true;
  readInfo.readSize    = 150;
  readInfo.interval    = GenomeInterval(0, 100, 250);
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

// Test the read fragment-size probability if that fragment supports SV.
// FragProb = min(C, 1-C),
// C = fragDist.cdf(read1Distance + read2Distance), where
//            fragDist = sample fragment distribution,
//            cdf = cumulative density function
//            read1Distance = distance from center pos of BP1 to read1
//            read2Distance = distance from center pos of BP2 to read2
// Following cases need to be tested:
// 1. If a fragment supports SV, frag probability is calculated as mentioned above
// 2. Fragment probability should be greater than 0.0001f.
// 3. Distance from breakpoint center position to read position should be greater than or equal to 50
// 4. For RNA FragProb = max(FragProb, 0.0001f) that means if FragProb is less than 0.0001f, then for RNA
//    always 0.0001f will be taken as fragment probability.
// 5. If read pair is inter chromosomal but breakpoint pair is not inter chromosomal, it throws an exception.
// 6. If Read1 chromosome id and bp1 chromosome id is not same, it throws an exception
// 7. If Read2 chromosome id and bp2 chromosome id is not same, it throws an exception
// 8. If Read1 is in reverse strand and sv.bp1.state == SVBreakendState::RIGHT_OPEN, it throws an exception
// 9. If Read2 is in reverse strand and sv.bp1.state == SVBreakendState::RIGHT_OPEN, it throws an exception
BOOST_AUTO_TEST_CASE(test_getFragProb)
{
  static const float eps(0.00000001);
  PairOptions        options1(false);  // false means it is a DNA sample
  bool               isFragSupportSV;
  float              fragProb;
  SVCandidate        candidate;
  // bp1 center pos = 254
  candidate.bp1.interval = GenomeInterval(0, 250, 260);
  candidate.bp1.state    = SVBreakendState::RIGHT_OPEN;
  // bp2 center pos = 276
  candidate.bp2.interval = GenomeInterval(0, 270, 280);
  candidate.bp2.state    = SVBreakendState::LEFT_OPEN;
  SVCandidateSetSequenceFragment fragment;
  // creating fragment distribution
  SizeDistribution fragDistribution;
  for (unsigned i(0); i < 250; ++i) {
    fragDistribution.addObservation(50);
    fragDistribution.addObservation(75);
    fragDistribution.addObservation(100);
    fragDistribution.addObservation(125);
  }
  std::string queryseq1 = "AGCTGACTGATCGATTTTTTACGTAGAGGAGCTTTGACGTATGAGCCTGATATGAGCCTG";
  std::string queryseq2 = "TGACGTATGAGCCTGATATGAGCCT";

  // This is mate-1 read which is in forward strand
  bam_record bamRecord1;
  buildTestBamRecord(bamRecord1, 0, 200, 0, 300, 60, 15, "60M", queryseq1, 125);
  bamRecord1.set_qname("Read-1");
  fragment.read1.bamrec = bamRecord1;

  // This is mate-2 read which is in reverse strand
  bam_record bamRecord2;
  buildTestBamRecord(bamRecord2, 0, 300, 0, 200, 25, 15, "25M", queryseq2, 125);
  bamRecord2.toggle_is_fwd_strand();
  bamRecord2.toggle_is_mate_fwd_strand();
  bamRecord2.set_qname("Read-2");
  fragment.read2.bamrec = bamRecord2;

  float expected(0.24999994);
  // Case-1, Case-2 and case-3 are desiged here.
  // Fragment-size probability = ~0.25 > 0.0001f
  // Distance between read-1 and BP1 center is 254 - 200 + 1 = 56 which is greater than 50
  // Distance between read-2 and BP2 center is 325 - 276 + 1 = 50 which is equal to 50.
  getFragProb(options1, candidate, fragment, fragDistribution, true, isFragSupportSV, fragProb);
  BOOST_REQUIRE(isFragSupportSV);
  BOOST_REQUIRE_CLOSE(fragProb, expected, eps);

  // This is similar test case as above but here fragment probality is performed
  // on RNA sample. So fragment-size probability is max(~0.25, 0.0001) = 0.25.
  PairOptions options2(true);
  getFragProb(options2, candidate, fragment, fragDistribution, true, isFragSupportSV, fragProb);
  BOOST_REQUIRE(isFragSupportSV);
  BOOST_REQUIRE_CLOSE(fragProb, expected, eps);

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
  // So Read-1 distance = 224 - 200 = 24 which is less than 50
  getFragProb(options1, candidate, fragment, fragDistribution, true, isFragSupportSV, fragProb);
  BOOST_REQUIRE(!isFragSupportSV);
  BOOST_REQUIRE_EQUAL(fragProb, 0);

  // Inter-chromosomal read pair
  bam_record bamRecord3;
  buildTestBamRecord(bamRecord3, 0, 200, 1, 300, 60, 15, "60M", queryseq1, 0);
  bamRecord3.set_qname("Read-3");
  fragment.read1.bamrec = bamRecord3;
  bam_record bamRecord4;
  buildTestBamRecord(bamRecord4, 1, 300, 0, 200, 25, 15, "25M", queryseq2, 0);
  bamRecord4.toggle_is_fwd_strand();
  bamRecord4.toggle_is_mate_fwd_strand();
  fragment.read2.bamrec = bamRecord4;
  // Breakpoint pair is not interchromosomal
  candidate.bp2.interval = GenomeInterval(0, 270, 280);
  // Case-5 are designed here where read pair is inter chromosomal
  // but breakpoint pair is not inter chromosomal
  BOOST_CHECK_THROW(
      getFragProb(options1, candidate, fragment, fragDistribution, true, isFragSupportSV, fragProb),
      illumina::common::GeneralException);

  // Read pair and breakpoint pair are inter chromosomal, but
  // chromosome of read-1 and BP-1 are not matching. Case-6 is designed here.
  candidate.bp1.interval = GenomeInterval(2, 250, 260);
  candidate.bp2.interval = GenomeInterval(0, 270, 280);
  BOOST_CHECK_THROW(
      getFragProb(options1, candidate, fragment, fragDistribution, true, isFragSupportSV, fragProb),
      illumina::common::GeneralException);

  // Read pair and breakpoint pair are inter chromosomal, but
  // chromosome of read-2 and BP-2 are not matching. Case-7 is designed here.
  candidate.bp1.interval = GenomeInterval(0, 250, 260);
  candidate.bp2.interval = GenomeInterval(2, 270, 280);
  BOOST_CHECK_THROW(
      getFragProb(options1, candidate, fragment, fragDistribution, true, isFragSupportSV, fragProb),
      illumina::common::GeneralException);

  // Case-8 is designed here where read orientation and breakpoint orientation are not matching
  bamRecord3.toggle_is_fwd_strand();
  bamRecord3.get_data()->core.mtid = 0;
  fragment.read1.bamrec            = bamRecord3;
  bamRecord4.set_target_id(0);
  fragment.read2.bamrec  = bamRecord4;
  candidate.bp1.interval = GenomeInterval(0, 250, 260);
  candidate.bp1.state    = SVBreakendState::RIGHT_OPEN;
  candidate.bp2.interval = GenomeInterval(0, 270, 280);
  candidate.bp2.state    = SVBreakendState::RIGHT_OPEN;
  // read1(bamRecord3) is in reverse strand and sv.bp1.state = SVBreakendState::RIGHT_OPEN
  BOOST_CHECK_THROW(
      getFragProb(options1, candidate, fragment, fragDistribution, true, isFragSupportSV, fragProb),
      illumina::common::GeneralException);

  // Case-9 is designed here where read orientation and breakpoint orientation are not matching
  // read1 is in forward strand
  bamRecord3.toggle_is_fwd_strand();
  fragment.read1.bamrec = bamRecord3;
  // read2(bamRecord4) is in reverse strand and sv.bp2.state == SVBreakendState::RIGHT_OPEN
  BOOST_CHECK_THROW(
      getFragProb(options1, candidate, fragment, fragDistribution, true, isFragSupportSV, fragProb),
      illumina::common::GeneralException);
}

// For each bam test whether a bam read satisfies the following criteria
// 1. Start of bam read should overlap with the search range where search range is calculated as:
//            search range start = BP center pos - (max fragment size of the sample - min fragment support)
//            search range end = BP center pos + (max fragment size of the sample - min fragment support) + 1
// 2. Fragment length should be greater than or equal to min fragment length of the sample
// 3. Fragment length should be less than or equal to max fragment length of the sample
// 4. minimum of (breakpoint center pos - fragment start + 1) and (fragment end - breakpoint center pos)
//    should be greater than minimum fragment threshold which is 50
// If the above 4 conditions are satisfied for a fragment, then that fragment supports allele on breakpoint.
// Here the test cases are constructed with bams.
BOOST_AUTO_TEST_CASE(test_processBamProcList)
{
  // whether alignment file tumor or normal
  const std::vector<bool>         bamFileInfo = {false, false};
  const bam_header_info           bamHeader(buildTestBamHeader());
  std::unique_ptr<SVLocusScanner> scanner(buildSVLocusScanner(bamHeader));
  const PairOptions               options1(false);
  SVCandidate                     candidate;
  candidate.insertSeq =
      "GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGT"
      "ATTTTCGTCTGGGGGGTGTGCACGCGATAGCATTGCGAGACGCTGGA";
  // BP1 center pos = 159
  candidate.bp1.interval.range = known_pos_range2(100, 220);
  // BP2 center pos = 309
  candidate.bp2.interval.range = known_pos_range2(250, 370);
  SVEvidence evidence;
  evidence.samples.resize(2);
  // As we are interested in BP1,
  // search range start = 159 - (125-50) = 84
  // search range end = 159 + (125-50) + 1 = 235
  // Test cases for search range are mentioned in test_nextBAMIndex in SVScorePairProcessorTest.cpp
  std::shared_ptr<SVScorePairRefProcessor> processor(
      new SVScorePairRefProcessor(bamFileInfo, scanner.operator*(), options1, candidate, true, evidence));
  std::vector<SVScorer::pairProcPtr> pairProcList = {processor};
  std::vector<SVScorer::streamPtr>   bamStreams   = {bamStream1, bamStream2};
  SVId                               id;
  SVEvidenceWriterData               svEvidenceWriterData(2);
  processBamProcList(bamStreams, id, pairProcList, svEvidenceWriterData);
  // Check info for 1st bam
  // bam read start = 9. It is not overlapping with search range [84, 235).
  // As a result of this, the fragment is not supporting allele on BP1.
  BOOST_REQUIRE(!evidence.getSampleEvidence(0)["bamRecord1"].ref.bp1.isFragmentSupport);
  // Bam read start = 109. It is overlapping with the search range. But
  // its fragment length is 49 which is less than minimum fragment length(50) of the sample
  // As a result of this, the fragment is not supporting allele on BP1.
  BOOST_REQUIRE(!evidence.getSampleEvidence(0)["bamRecord2"].ref.bp1.isFragmentSupport);
  // Bam read start = 109. It is overlapping with the search range. But
  // its fragment length is 130 which is greater than maximum fragment length(125) of the
  // sample. As a result of this, the fragment is not supporting allele on BP1.
  BOOST_REQUIRE(!evidence.getSampleEvidence(0)["bamRecord3"].ref.bp1.isFragmentSupport);
  // Here min(159-109+1, 208-159) = 51 which is greater than 50
  // And also case-1, case-2 and case-3 are satisfied.
  // Fragment start = 109 which is overlapping with the search range[84,235).
  // Fragment size = 100 which is greater than 50 and less than 125.
  // As a result of this, fragment support of this
  // bam read will be set.
  BOOST_REQUIRE(evidence.getSampleEvidence(0)["bamRecord4"].ref.bp1.isFragmentSupport);

  // Check info for 2nd bam
  // bam read start = 10. It is not overlapping with search range [84, 235).
  // As a result of this, the fragment is not supporting allele on BP1.
  BOOST_REQUIRE(!evidence.getSampleEvidence(1)["bamRecord5"].ref.bp1.isFragmentSupport);
  // Bam read start = 108. It does overlap with the search range. However,
  // its fragment length is 49 which is less than minimum fragment length(50) of the sample
  // As a result of this, the fragment is not supporting allele on BP1.
  BOOST_REQUIRE(!evidence.getSampleEvidence(1)["bamRecord6"].ref.bp1.isFragmentSupport);
  // Bam read start = 108. It does overlap with the search range. However,
  // its fragment length is 130 which is greater than maximum fragment length(125) of the
  // sample. As a result of this, the fragment is not supporting allele on BP1.
  BOOST_REQUIRE(!evidence.getSampleEvidence(1)["bamRecord7"].ref.bp1.isFragmentSupport);
  // Here min(159-109+1, 208-159) = 51 which is greater than 50
  // And also case-1, case-2 and case-3 are satisfied.
  // Fragment start = 109 which is overlapping with the search range[84,235).
  // Fragment size = 100 which is greater than 50 and less than 125.
  // As a result of this, fragment support of this
  // bam read will be set.
  BOOST_REQUIRE(evidence.getSampleEvidence(1)["bamRecord8"].ref.bp1.isFragmentSupport);
}

// For Each sample Api computes the following informations:
// 1) Whether read-1 and read-2 are scanned
// 2) Whether read-1 and read-2 are anchored reads
// 3) Calculate the fragment-size probability as:
//       FragProb = min(C, 1-C),
//       C = fragDist.cdf(read1Distance + read2Distance), where
//            fragDist = sample fragment distribution,
//            cdf = cumulative density function
//            read1Distance = distance from center pos of BP1 to read1
//            read2Distance = distance from center pos of BP2 to read2
// 4) If the fragment supports allele on the BPs, then set that flag.
BOOST_AUTO_TEST_CASE(test_processExistingAltPairInfo)
{
  static const float              eps(0.00000001);
  const bam_header_info           bamHeader(buildTestBamHeader());
  std::unique_ptr<SVLocusScanner> scanner(buildSVLocusScanner(bamHeader));
  GSCOptions                      options;
  options.alignFileOpt.alignmentFilenames = {bamFileName1, bamFileName2};
  options.alignFileOpt.isAlignmentTumor   = {false, false};
  SVScorer     scorer(options, scanner.operator*(), bamHeader);
  TestSVScorer fSVScorer;
  PairOptions  pairOptions(false);  // false means it is a DNA sample

  // Candidate SV
  SVCandidate candidate;
  candidate.bp1.interval = GenomeInterval(0, 250, 260);
  candidate.bp1.state    = SVBreakendState::RIGHT_OPEN;
  candidate.bp2.interval = GenomeInterval(0, 270, 280);
  candidate.bp2.state    = SVBreakendState::LEFT_OPEN;

  // Creating a read pair for 1st bam file
  std::string queryseq1 = "AGCTGACTGATCGATTTTTTACGTAGAGGAGCTTTGACGTATGAGCCTGATATGAGCCTG";
  std::string queryseq2 = "TGACGTATGAGCCTGATATGAGCCT";
  bam_record  bamRecord1;
  buildTestBamRecord(bamRecord1, 0, 200, 0, 300, 60, 20, "60M", queryseq1, 125);
  bamRecord1.set_qname("Read-1");
  bamRecord1.toggle_is_first();
  bam_record bamRecord2;
  buildTestBamRecord(bamRecord2, 0, 300, 0, 200, 25, 10, "25M", queryseq2, 125);
  bamRecord2.toggle_is_second();
  bamRecord2.toggle_is_fwd_strand();
  bamRecord2.toggle_is_mate_fwd_strand();
  bamRecord2.set_qname("Read-1");

  // creating another read pair for 2nd bam file
  bam_record bamRecord3;
  buildTestBamRecord(bamRecord3, 0, 200, 0, 310, 60, 10, "60M", queryseq1, 125);
  bamRecord3.set_qname("Read-2");
  bamRecord3.toggle_is_first();
  bam_record bamRecord4;
  buildTestBamRecord(bamRecord4, 0, 310, 0, 200, 25, 20, "25M", queryseq2, 125);
  bamRecord4.toggle_is_second();
  bamRecord4.toggle_is_fwd_strand();
  bamRecord4.toggle_is_mate_fwd_strand();
  bamRecord4.set_qname("Read-2");

  // Creating fragment informations for 1st bam file
  SVCandidateSetData                         candidateSetData;
  SVCandidateSetSequenceFragmentSampleGroup& group1 = candidateSetData.getDataGroup(0);
  group1.add(bamHeader, bamRecord1, false, true, true);
  group1.add(bamHeader, bamRecord2, false, true, true);
  SVSequenceFragmentAssociation association(0, SVEvidenceType::PAIR);
  group1.begin()->svLink.push_back(association);
  // Creating fragment information for 2nd bam
  SVCandidateSetSequenceFragmentSampleGroup& group2 = candidateSetData.getDataGroup(1);
  group2.add(bamHeader, bamRecord3, false, true, true);
  group2.add(bamHeader, bamRecord4, false, true, true);
  group2.begin()->svLink.push_back(association);

  SVId       id;
  SVEvidence evidence;
  evidence.samples.resize(2);
  SVEvidenceWriterData svEvidenceWriterData(2);
  fSVScorer.processExistingAltPairInfo(
      scorer, pairOptions, candidateSetData, candidate, id, evidence, svEvidenceWriterData);

  float expected(0.24999994);
  // Check 1st bam file
  for (const SVCandidateSetSequenceFragment& fragment : group1) {
    SVFragmentEvidence       fragmentEvidence(evidence.getSampleEvidence(0)[fragment.qname()]);
    SVFragmentEvidenceAllele alt(fragmentEvidence.alt);
    // Whether read-1 is scanned
    BOOST_REQUIRE(fragmentEvidence.read1.isScanned);
    // Whether read-2 is scanned
    BOOST_REQUIRE(fragmentEvidence.read2.isScanned);
    // bamRecord-1 is anchored read as mapping quality (20) is
    // greater than min mappig quality (15)
    BOOST_REQUIRE(fragmentEvidence.read1.isAnchored(false));
    // bamRecord-2 is not anchored read as mapping quality (10)
    // is less than min mappig quality(15)
    BOOST_REQUIRE(!fragmentEvidence.read2.isAnchored(false));
    // FragProb = min(C, 1-C),
    // C = fragDist.cdf(read1Distance + read2Distance), where
    //            fragDist = sample fragment distribution,
    //            cdf = cumulative density function
    //            read1Distance = distance from center pos of BP1 to read1
    //            read2Distance = distance from center pos of BP2 to read2
    // Detail has been explained in test_getFragProb
    BOOST_REQUIRE_CLOSE(alt.bp1.fragLengthProb, expected, eps);
    BOOST_REQUIRE_CLOSE(alt.bp2.fragLengthProb, expected, eps);
    // Fragment is supporting alt allele in BP1 and BP2
    BOOST_REQUIRE(alt.bp1.isFragmentSupport);
    BOOST_REQUIRE(alt.bp2.isFragmentSupport);
  }

  // Check 2nd bam file
  for (const SVCandidateSetSequenceFragment& fragment : group2) {
    SVFragmentEvidence       fragmentEvidence(evidence.getSampleEvidence(1)[fragment.qname()]);
    SVFragmentEvidenceAllele alt(fragmentEvidence.alt);
    // Whether read-1 is scanned
    BOOST_REQUIRE(fragmentEvidence.read1.isScanned);
    // Whether read-2 is scanned
    BOOST_REQUIRE(fragmentEvidence.read2.isScanned);
    // bamRecord-3 is not anchored read as mapping quality (10) is
    // less than min mappig quality (15)
    BOOST_REQUIRE(!fragmentEvidence.read1.isAnchored(false));
    // bamRecord-4 is anchored read as mapping quality (20) is
    // greater than min mappig quality (15)
    BOOST_REQUIRE(fragmentEvidence.read2.isAnchored(false));
    // FragProb = min(C, 1-C),
    // C = fragDist.cdf(read1Distance + read2Distance), where
    //            fragDist = sample fragment distribution,
    //            cdf = cumulative density function
    //            read1Distance = distance from center pos of BP1 to read1
    //            read2Distance = distance from center pos of BP2 to read2
    // Detail has been explained in test_getFragProb
    BOOST_REQUIRE_CLOSE(alt.bp1.fragLengthProb, expected, eps);
    BOOST_REQUIRE_CLOSE(alt.bp2.fragLengthProb, expected, eps);
    // Fragment is supporting alt allele in BP1 and BP2
    BOOST_REQUIRE(alt.bp1.isFragmentSupport);
    BOOST_REQUIRE(alt.bp2.isFragmentSupport);
  }
}

// Following two cases are designed here:
// Let's say, 95th percentile value is F and
// deletion size from candidate SV is D.
// 1. If D <= 2*F then it is incomplete alt pair info, so it will not support allele on BPs
// 2. If D > 2*F, then fragment probability is calculated as:
//       FragProb = min(C, 1-C),
//       C = fragDist.cdf(read1Distance + read2Distance), where
//            fragDist = sample fragment distribution,
//            cdf = cumulative density function
//            read1Distance = distance from center pos of BP1 to read1
//            read2Distance = distance from center pos of BP2 to read2
BOOST_AUTO_TEST_CASE(test_getSVPairSupport)
{
  static const float eps(0.00000001);
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
  jumpAlignmentResultType.align1    = alignment1;
  jumpAlignmentResultType.align2    = alignment2;
  jumpAlignmentResultType.jumpRange = 2;
  candidateAssemblyData.spanningAlignments.push_back(jumpAlignmentResultType);
  candidateAssemblyData.isSpanning = true;
  candidateAssemblyData.extendedContigs.push_back(
      "GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGT"
      "ATTTTCGTCTGGGGGGTGTGCACGCGATAGCATTGCGAGACGCTGGA");
  candidateAssemblyData.isCandidateSpanning = true;
  // Test SV  locus created where 95th percentile value is 125
  const bam_header_info           bamHeader(buildTestBamHeader());
  std::unique_ptr<SVLocusScanner> scanner(buildTestSVLocusScanner(bamHeader));

  GSCOptions options;
  options.alignFileOpt.alignmentFilenames = {bamFileName1};
  options.alignFileOpt.isAlignmentTumor   = {false};
  SVScorer     scorer(options, scanner.operator*(), bamHeader);
  TestSVScorer fSVScorer;
  SVCandidate  candidate1;
  candidate1.bp1.interval = GenomeInterval(0, 250, 260);
  candidate1.bp1.state    = SVBreakendState::RIGHT_OPEN;
  candidate1.bp2.interval = GenomeInterval(0, 270, 280);
  candidate1.bp2.state    = SVBreakendState::LEFT_OPEN;
  std::string queryseq1   = "AGCTGACTGATCGATTTTTTACGTAGAGGAGCTTTGACGTATGAGCCTGATATGAGCCTG";
  std::string queryseq2   = "TGACGTATGAGCCTGATATGAGCCT";
  bam_record  bamRecord1;
  buildTestBamRecord(bamRecord1, 0, 200, 0, 300, 60, 15, "60M", queryseq1, 125);
  bamRecord1.set_qname("Read-1");
  bamRecord1.toggle_is_first();
  bam_record bamRecord2;
  buildTestBamRecord(bamRecord2, 0, 300, 0, 200, 25, 15, "25M", queryseq2, 125);
  bamRecord2.toggle_is_second();
  bamRecord2.toggle_is_fwd_strand();
  bamRecord2.toggle_is_mate_fwd_strand();
  bamRecord2.set_qname("Read-1");

  // Creating fragment informations
  SVCandidateSetData                         candidateSetData1;
  SVCandidateSetSequenceFragmentSampleGroup& group1 = candidateSetData1.getDataGroup(0);
  group1.add(bamHeader, bamRecord1, false, true, true);
  group1.add(bamHeader, bamRecord2, false, true, true);
  SVSequenceFragmentAssociation association(0, SVEvidenceType::PAIR);
  group1.begin()->svLink.push_back(association);
  SVId       id;
  SVEvidence evidence1;
  evidence1.samples.resize(1);
  SVEvidenceWriterData svEvidenceWriterData1(1);
  candidate1.setPrecise();
  candidate1.assemblyAlignIndex = 0;
  fSVScorer.getSVPairSupport(
      scorer, candidateSetData1, candidateAssemblyData, candidate1, id, evidence1, svEvidenceWriterData1);

  // According to description, F = 125 and D = 270-250 = 20
  // So D < 2*F (20 < 250). So fragment is not supporting allele on the BPs.
  // For that reason fragment-size probability is zero.
  for (const SVCandidateSetSequenceFragment& fragment : group1) {
    SVFragmentEvidenceAllele alt(evidence1.getSampleEvidence(0)[fragment.qname()].alt);
    BOOST_REQUIRE_EQUAL(alt.bp1.fragLengthProb, 0);
    BOOST_REQUIRE_EQUAL(alt.bp2.fragLengthProb, 0);
    BOOST_REQUIRE(!alt.bp1.isFragmentSupport);
    BOOST_REQUIRE(!alt.bp2.isFragmentSupport);
  }

  SVCandidate candidate2;
  candidate2.bp1.interval = GenomeInterval(0, 250, 260);
  candidate2.bp1.state    = SVBreakendState::RIGHT_OPEN;
  candidate2.bp2.interval = GenomeInterval(0, 770, 780);
  candidate2.bp2.state    = SVBreakendState::LEFT_OPEN;
  bam_record bamRecord3;
  buildTestBamRecord(bamRecord3, 0, 200, 0, 800, 60, 15, "60M", queryseq1, 625);
  bamRecord3.set_qname("Read-2");
  bamRecord3.toggle_is_first();

  bam_record bamRecord4;
  buildTestBamRecord(bamRecord4, 0, 800, 0, 200, 25, 15, "25M", queryseq2, 625);
  bamRecord4.toggle_is_second();
  bamRecord4.toggle_is_fwd_strand();
  bamRecord4.toggle_is_mate_fwd_strand();
  bamRecord4.set_qname("Read-2");

  SVCandidateSetData                         candidateSetData2;
  SVCandidateSetSequenceFragmentSampleGroup& group2 = candidateSetData2.getDataGroup(0);
  group2.add(bamHeader, bamRecord3, false, true, true);
  group2.add(bamHeader, bamRecord4, false, true, true);
  group2.begin()->svLink.push_back(association);
  SVEvidence evidence2;
  evidence2.samples.resize(1);
  SVEvidenceWriterData svEvidenceWriterData2(1);
  candidate2.setPrecise();
  candidate2.assemblyAlignIndex = 0;
  float expected                = 0.24999994;
  fSVScorer.getSVPairSupport(
      scorer, candidateSetData2, candidateAssemblyData, candidate2, id, evidence2, svEvidenceWriterData2);
  // Here F = 125, But D = 770-250 = 520 which is greater than 2*F(250).
  // So fragment is supporting alt allele on the BPs.
  for (const SVCandidateSetSequenceFragment& fragment : group2) {
    SVFragmentEvidenceAllele alt(evidence2.getSampleEvidence(0)[fragment.qname()].alt);
    BOOST_REQUIRE_CLOSE(alt.bp1.fragLengthProb, expected, eps);
    BOOST_REQUIRE_CLOSE(alt.bp2.fragLengthProb, expected, eps);
    BOOST_REQUIRE(alt.bp1.isFragmentSupport);
    BOOST_REQUIRE(alt.bp2.isFragmentSupport);
  }
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
