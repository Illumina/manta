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

#include "boost/test/unit_test.hpp"

#include "manta/SVLocusEvidenceCount.hpp"
#include "test/testAlignmentDataUtil.hpp"
#include "test/testFileMakers.hpp"
#include "test/testSVLocusScanner.hpp"
#include "test/testSVLocusUtil.hpp"
#include "test/testUtil.hpp"

#include "SVFinder.cpp"
#include "SVFinder.hpp"

BOOST_AUTO_TEST_SUITE(SVFinderTest_test_suite)

// Test the fraction of anomalous or split evidence count to total evidence count
BOOST_AUTO_TEST_CASE(test_SpanningNoiseRate)
{
  AllSampleReadCounts counts;
  counts.setSampleCount(2);
  SampleReadCounts sample1(counts.getSampleCounts(0));
  sample1.input.evidenceCount.anom         = 10;
  sample1.input.evidenceCount.split        = 5;
  sample1.input.evidenceCount.anomAndSplit = 4;
  sample1.input.evidenceCount.total        = 19;
  counts.getSampleCounts(0).merge(sample1);

  SampleReadCounts sample2(counts.getSampleCounts(1));
  sample2.input.evidenceCount.anom         = 25;
  sample2.input.evidenceCount.split        = 5;
  sample2.input.evidenceCount.anomAndSplit = 10;
  sample2.input.evidenceCount.total        = 40;
  counts.getSampleCounts(1).merge(sample2);

  BOOST_REQUIRE_EQUAL(getSpanningNoiseRate(counts, 0), 0.020608439646712464);
  BOOST_REQUIRE_EQUAL(getSpanningNoiseRate(counts, 1), 0.028846153846153848);
}

// Test the fraction of semi-aligned evidence count to total evidence count
BOOST_AUTO_TEST_CASE(test_AssemblyNoiseRate)
{
  AllSampleReadCounts counts;
  counts.setSampleCount(2);
  SampleReadCounts sample1(counts.getSampleCounts(0));
  sample1.input.evidenceCount.assm  = 10;
  sample1.input.evidenceCount.total = 19;
  counts.getSampleCounts(0).merge(sample1);

  SampleReadCounts sample2(counts.getSampleCounts(1));
  sample2.input.evidenceCount.assm  = 25;
  sample2.input.evidenceCount.total = 40;
  counts.getSampleCounts(1).merge(sample2);

  BOOST_REQUIRE_EQUAL(getAssemblyNoiseRate(counts, 0), 0.019627085377821395);
  BOOST_REQUIRE_EQUAL(getAssemblyNoiseRate(counts, 1), 0.033653846153846152);
}

// test if read supports an SV on this edge, if so, add to SVData
BOOST_AUTO_TEST_CASE(test_AddSVNodeRead)
{
  const bam_header_info           bamHeader(buildTestBamHeader());
  std::unique_ptr<SVLocusScanner> scanner(buildTestSVLocusScanner(bamHeader));
  SampleEvidenceCounts            eCounts;

  const unsigned                 defaultReadGroupIndex(0);
  const reference_contig_segment refSeq;

  // supplementary read in SV evidence
  bam_record supplementSASplitRead;
  buildTestBamRecord(supplementSASplitRead);                 // pos=100, matePos=200 and fragmentSize=100
  addSupplementaryAlignmentEvidence(supplementSASplitRead);  // SA tag : chrFoo,300,-,54H22M,50,0;

  SVLocus locus1;
  locus1.addNode(GenomeInterval(0, 80, 120));
  locus1.addNode(GenomeInterval(0, 279, 319));
  locus1.addNode(GenomeInterval(0, 410, 450));

  SVCandidateSetSequenceFragmentSampleGroup svDatagroup;

  // test a supplementSASplitRead read is overlapping with a locus node when localnode's coordinate is
  // GenomeInterval(0,80,120) and remoteNode's coordinate is GenomeInterval(0,279,319). It will add an entry
  // in svDatagroup. As supplementSASplitRead is a split read, so the evidence genomic intervals found from
  // supplementSASplitRead are
  // as follows:
  // 1. Local Evidence = GenomeInterval: 0:[80,120) (Read posotion is 100 and 20 bp padding on each side)
  // 2. Remote Evidence = GenomeInterval: 0:[279,319) (According to SA tag other segment position is 299 and
  // 20bp
  //                                                   padding on each side)
  addSVNodeRead(
      bamHeader,
      scanner.operator*(),
      ((const SVLocus&)locus1).getNode(0),
      ((const SVLocus&)locus1).getNode(1),
      supplementSASplitRead,
      defaultReadGroupIndex,
      true,
      refSeq,
      true,
      false,
      svDatagroup,
      eCounts);
  BOOST_REQUIRE_EQUAL(svDatagroup.size(), 1u);

  // This read is part of the anomalous read pair although there is a large insertion,
  // priority always goes to anomalous read pair. Based on the read's start-end and its
  // mate's start-end, the evidence genomic intervals are as follows :
  // 1. Local Evidence = GenomeInterval: 0:[400,440) (Based on read end and 40 bp minPairBreakendSize)
  // 2. Remote Evidence = GenomeInterval: 0:[260,300) (Based on mate's start and 40 bp minPairBreakendSize)
  bam_record anomalousReadBasedOnFragment;
  buildTestBamRecord(anomalousReadBasedOnFragment, 0, 200, 0, 300, 2200, 15, "100M2000I100M");
  anomalousReadBasedOnFragment.set_qname("anomalous_read");

  // test a read is not overlapping with a locus node (mentioned above) when localnode's coordinate is
  // GenomeInterval(0,80,120) and remoteNode's coordinate is GenomeInterval(0,279,319). It will not add any
  // entry in svDatagroup.
  addSVNodeRead(
      bamHeader,
      scanner.operator*(),
      ((const SVLocus&)locus1).getNode(0),
      ((const SVLocus&)locus1).getNode(1),
      anomalousReadBasedOnFragment,
      defaultReadGroupIndex,
      true,
      refSeq,
      true,
      false,
      svDatagroup,
      eCounts);
  BOOST_REQUIRE_EQUAL(svDatagroup.size(), 1u);

  // test a read is overlapping with the locus node (mentioned above) when localnode's coordinate is
  // GenomeInterval(0,410,450) and remoteNode's coordinate is GenomeInterval(0,279,319). It will add another
  // entry in svDatagroup.
  addSVNodeRead(
      bamHeader,
      scanner.operator*(),
      ((const SVLocus&)locus1).getNode(2),
      ((const SVLocus&)locus1).getNode(1),
      anomalousReadBasedOnFragment,
      defaultReadGroupIndex,
      true,
      refSeq,
      true,
      false,
      svDatagroup,
      eCounts);
  BOOST_REQUIRE_EQUAL(svDatagroup.size(), 2u);
}

// test reference sequence of a segment. It will add 100 bases on both side
// that means if Genomic start and end coordinates are 1 and 2 respectively, and chromosome id
// is 0,  then the modified interval will be [max(0, 1-100), min(2+100, chrLength)).
// So the total length will be 102.
BOOST_AUTO_TEST_CASE(test_GetNodeRef)
{
  const bam_header_info bamHeader(buildTestBamHeader());
  SVLocus               locus;
  locus.addNode(GenomeInterval(0, 1, 2));
  GenomeInterval           searchInterval;
  reference_contig_segment refSeq;
  getNodeRefSeq(bamHeader, locus, 0, getTestReferenceFilename(), searchInterval, refSeq);
  // check the size first
  BOOST_REQUIRE_EQUAL(refSeq.seq().size(), 102);
  // check the sequence
  BOOST_REQUIRE_EQUAL(
      refSeq.seq(),
      "GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGG"
      "TATTTTCGTCTGGGGGGTGTGCACGCGATAGCATTGCGAGACGCTGGA");
}

// test candidates must have at least evidence of 2
BOOST_AUTO_TEST_CASE(test_IsCandidateCountSufficient)
{
  SVCandidate candidate;
  for (unsigned i(0); i < SVEvidenceType::SIZE; ++i) candidate.bp1.lowresEvidence.add(i, 1);

  // Evidence count is not sufficient
  BOOST_REQUIRE(!isCandidateCountSufficient(candidate));

  for (unsigned i(0); i < SVEvidenceType::SIZE; ++i) candidate.bp1.lowresEvidence.add(i, 1);

  // Evidence count is sufficient
  BOOST_REQUIRE(isCandidateCountSufficient(candidate));
}

// test depth on each location i.e. number of
// read bases overlap in a location.
BOOST_AUTO_TEST_CASE(test_AddReadToDepthEst)
{
  bam_record bamRecord1;
  buildTestBamRecord(bamRecord1, 0, 200, 0, 210, 20, 15);
  bamRecord1.set_qname("Read-1");

  bam_record bamRecord2;
  buildTestBamRecord(bamRecord2, 0, 210, 0, 220, 20, 15);
  bamRecord2.set_qname("Read-2");

  std::vector<unsigned> depth(30);
  addReadToDepthEst(bamRecord1, 200, depth);
  addReadToDepthEst(bamRecord2, 200, depth);

  // test the coverage
  for (unsigned i = 0; i < 30; i++) {
    if (i >= 10 && i <= 19)
      BOOST_REQUIRE_EQUAL(depth[i], 2u);  // second bamRead starts 10 bases after first bamRead
    else
      BOOST_REQUIRE_EQUAL(depth[i], 1u);
  }
}

// test the significance of a break point based on the supporting read
// observations relative to a background noise rate.
BOOST_AUTO_TEST_CASE(test_IsBreakPointSignificant)
{
  std::vector<double> signalReadInfo;

  // minimum signal count should be 2
  BOOST_REQUIRE(!isBreakPointSignificant(0.1, 0.5, signalReadInfo));

  signalReadInfo.push_back(96);
  signalReadInfo.push_back(158);
  signalReadInfo.push_back(163);
  // Break point is not significant.
  // 0.005 is alpha in the binomial significance test.
  BOOST_REQUIRE(!isBreakPointSignificant(0.005, 0.005, signalReadInfo));

  signalReadInfo.clear();
  signalReadInfo.push_back(3440);
  signalReadInfo.push_back(3443);
  signalReadInfo.push_back(3452);
  signalReadInfo.push_back(3489);
  // Break point is significant
  // 0.03 is alpha in the binomial significance test.
  BOOST_REQUIRE(isBreakPointSignificant(0.03, 0.008, signalReadInfo));
}

// test the significance of a spanning candidate for  minimum supporting evidence.
// Sppanning candidate is significant if either break point 1 or break point 2 is
// significant. This test verifies following cases:
// 1) When no breakpoint is significant.
// 2) When Breakpoint-1 is significant and Breakpoint-2 is not significant.
// 3) When Breakpoint-2 is significant and Breakpoint-1 is not significant.
// 4) When both the Breakpoints are significant.
BOOST_AUTO_TEST_CASE(test_IsSpanningCandidateSignalSignificant)
{
  SVCandidate    svCandidate;
  FatSVCandidate fatSVCandidate(svCandidate, 1u);
  // Spanning candidate is not significant as none of the breakpoint
  // satisfies minimum evidence(2) criteria.
  BOOST_REQUIRE(!isSpanningCandidateSignalSignificant(0.008, fatSVCandidate, 0));

  // test when both breakpoint-1 and breakpoint-2 are not significant where
  // 0.03 is alpha in the binomial signicance test.
  // spanning noise rate is 0.008 which is the probability of success in
  // binomial significance test.
  fatSVCandidate.bp1EvidenceIndex[0][0].push_back(3443);
  fatSVCandidate.bp1EvidenceIndex[0][0].push_back(3468);
  fatSVCandidate.bp1EvidenceIndex[0][0].push_back(3520);
  fatSVCandidate.bp1EvidenceIndex[0][0].push_back(3569);

  fatSVCandidate.bp2EvidenceIndex[0][0].push_back(1403);
  fatSVCandidate.bp2EvidenceIndex[0][0].push_back(1428);
  fatSVCandidate.bp2EvidenceIndex[0][0].push_back(1480);
  fatSVCandidate.bp2EvidenceIndex[0][0].push_back(1507);
  BOOST_REQUIRE(!isSpanningCandidateSignalSignificant(0.008, fatSVCandidate, 0));

  // test when breakpoint-1 is significant and
  // breakpoint-2 is not significant.
  fatSVCandidate.bp1EvidenceIndex[0][0].clear();
  fatSVCandidate.bp2EvidenceIndex[0][0].clear();
  fatSVCandidate.bp1EvidenceIndex[0][0].push_back(3443);
  fatSVCandidate.bp1EvidenceIndex[0][0].push_back(3452);
  fatSVCandidate.bp1EvidenceIndex[0][0].push_back(3440);
  fatSVCandidate.bp1EvidenceIndex[0][0].push_back(3489);

  fatSVCandidate.bp2EvidenceIndex[0][0].push_back(1403);
  fatSVCandidate.bp2EvidenceIndex[0][0].push_back(1428);
  fatSVCandidate.bp2EvidenceIndex[0][0].push_back(1480);
  fatSVCandidate.bp2EvidenceIndex[0][0].push_back(1507);
  BOOST_REQUIRE(isSpanningCandidateSignalSignificant(0.008, fatSVCandidate, 0));

  // test when breakpoint-2 is significant and
  // breakpoint-1 is not significant.
  fatSVCandidate.bp1EvidenceIndex[0][0].clear();
  fatSVCandidate.bp2EvidenceIndex[0][0].clear();
  fatSVCandidate.bp1EvidenceIndex[0][0].push_back(1403);
  fatSVCandidate.bp1EvidenceIndex[0][0].push_back(1428);
  fatSVCandidate.bp1EvidenceIndex[0][0].push_back(1480);
  fatSVCandidate.bp1EvidenceIndex[0][0].push_back(1507);

  fatSVCandidate.bp2EvidenceIndex[0][0].push_back(3443);
  fatSVCandidate.bp2EvidenceIndex[0][0].push_back(3452);
  fatSVCandidate.bp2EvidenceIndex[0][0].push_back(3440);
  fatSVCandidate.bp2EvidenceIndex[0][0].push_back(3489);
  BOOST_REQUIRE(isSpanningCandidateSignalSignificant(0.008, fatSVCandidate, 0));

  // test when both breakpoint-1 and breakpoint-2 are significant.
  fatSVCandidate.bp1EvidenceIndex[0][0].clear();
  fatSVCandidate.bp2EvidenceIndex[0][0].clear();
  fatSVCandidate.bp1EvidenceIndex[0][0].push_back(3443);
  fatSVCandidate.bp1EvidenceIndex[0][0].push_back(3452);
  fatSVCandidate.bp1EvidenceIndex[0][0].push_back(3440);
  fatSVCandidate.bp1EvidenceIndex[0][0].push_back(3489);

  fatSVCandidate.bp2EvidenceIndex[0][0].push_back(1403);
  fatSVCandidate.bp2EvidenceIndex[0][0].push_back(1412);
  fatSVCandidate.bp2EvidenceIndex[0][0].push_back(1400);
  fatSVCandidate.bp2EvidenceIndex[0][0].push_back(1449);
  BOOST_REQUIRE(isSpanningCandidateSignalSignificant(0.008, fatSVCandidate, 0));
}

// test the significance of a complex candidate for  minimum supporting evidence
// where complex means that we have no specific hypothesis for the SV -
// it is just a single genomic region for which we schedule local assembly.
BOOST_AUTO_TEST_CASE(test_IsComplexCandidateSignalSignificant)
{
  SVCandidate    svCandidate;
  FatSVCandidate fatSVCandidate(svCandidate, 1u);

  // Complex break point is not significant where
  // assembly noise rate is 0.008 which is the probability of success in
  // binomial significance test.
  fatSVCandidate.bp1EvidenceIndex[0][0].push_back(3443);
  fatSVCandidate.bp1EvidenceIndex[0][0].push_back(3452);
  fatSVCandidate.bp1EvidenceIndex[0][0].push_back(3440);
  fatSVCandidate.bp1EvidenceIndex[0][0].push_back(3489);
  BOOST_REQUIRE(!isComplexCandidateSignalSignificant(0.008, fatSVCandidate, 0));

  // Complex break point is significant
  fatSVCandidate.bp1EvidenceIndex[0][0].clear();
  fatSVCandidate.bp1EvidenceIndex[0][0].push_back(3443);
  fatSVCandidate.bp1EvidenceIndex[0][0].push_back(3452);
  fatSVCandidate.bp1EvidenceIndex[0][0].push_back(3440);
  fatSVCandidate.bp1EvidenceIndex[0][0].push_back(3448);
  BOOST_REQUIRE(isComplexCandidateSignalSignificant(0.008, fatSVCandidate, 0));
}

// test the significance of spanning candidate across all the bams
// relative to spanning noise rate. This test checks whether the method
// returns true if one of the bams has significant spanning candidate.
BOOST_AUTO_TEST_CASE(test_IsAnySpanningCandidateSignalSignificant)
{
  SVCandidate    svCandidate;
  FatSVCandidate fatSVCandidate(svCandidate, 2u);  // fat sv candidate object for 2 bams

  // insert read index values for 1st bam
  fatSVCandidate.bp1EvidenceIndex[0][0].push_back(3443);
  fatSVCandidate.bp1EvidenceIndex[0][0].push_back(3452);
  fatSVCandidate.bp1EvidenceIndex[0][0].push_back(3440);
  fatSVCandidate.bp1EvidenceIndex[0][0].push_back(3489);

  fatSVCandidate.bp2EvidenceIndex[0][0].push_back(1403);
  fatSVCandidate.bp2EvidenceIndex[0][0].push_back(1428);
  fatSVCandidate.bp2EvidenceIndex[0][0].push_back(1480);
  fatSVCandidate.bp2EvidenceIndex[0][0].push_back(1507);

  // insert read index values for 2nd bam
  fatSVCandidate.bp1EvidenceIndex[0][1].push_back(3443);
  fatSVCandidate.bp1EvidenceIndex[0][1].push_back(3452);
  fatSVCandidate.bp1EvidenceIndex[0][1].push_back(3440);
  fatSVCandidate.bp1EvidenceIndex[0][1].push_back(3489);

  fatSVCandidate.bp2EvidenceIndex[0][1].push_back(1403);
  fatSVCandidate.bp2EvidenceIndex[0][1].push_back(1428);
  fatSVCandidate.bp2EvidenceIndex[0][1].push_back(1480);
  fatSVCandidate.bp2EvidenceIndex[0][1].push_back(1507);

  std::vector<double> spanningNoiseRate;
  spanningNoiseRate.push_back(0.008);  // 1st bam spanning noise rate
  spanningNoiseRate.push_back(0.1);    // 2nd bam spanning noise rate

  // spanning candidate is significant for 1st bam
  BOOST_REQUIRE(isAnySpanningCandidateSignalSignificant(1, fatSVCandidate, spanningNoiseRate));

  spanningNoiseRate.clear();
  spanningNoiseRate.push_back(0.1);  // 1st bam spanning noise rate
  spanningNoiseRate.push_back(0.1);  // 2nd bam spanning noise rate
  // spanning candidate is not significant for any of the bams
  BOOST_REQUIRE(!isAnySpanningCandidateSignalSignificant(1, fatSVCandidate, spanningNoiseRate));
}

// test the significance of complex candidate across all the bams
// relative to assembly noise rate. This test checks whether the method
// returns true if one of the bams has complex candidate.
BOOST_AUTO_TEST_CASE(test_IsAnyComplexCandidateSignalSignificant)
{
  SVCandidate    svCandidate;
  FatSVCandidate fatSVCandidate(svCandidate, 2u);  // fat sv candidate object for 2 bams

  // insert read index values for 1st bam
  fatSVCandidate.bp1EvidenceIndex[0][0].push_back(3443);
  fatSVCandidate.bp1EvidenceIndex[0][0].push_back(3452);
  fatSVCandidate.bp1EvidenceIndex[0][0].push_back(3440);
  fatSVCandidate.bp1EvidenceIndex[0][0].push_back(3489);

  fatSVCandidate.bp2EvidenceIndex[0][0].push_back(1403);
  fatSVCandidate.bp2EvidenceIndex[0][0].push_back(1428);
  fatSVCandidate.bp2EvidenceIndex[0][0].push_back(1480);
  fatSVCandidate.bp2EvidenceIndex[0][0].push_back(1507);

  // insert values for 2nd bam
  fatSVCandidate.bp1EvidenceIndex[0][1].push_back(3443);
  fatSVCandidate.bp1EvidenceIndex[0][1].push_back(3452);
  fatSVCandidate.bp1EvidenceIndex[0][1].push_back(3440);
  fatSVCandidate.bp1EvidenceIndex[0][1].push_back(3489);

  fatSVCandidate.bp2EvidenceIndex[0][1].push_back(1403);
  fatSVCandidate.bp2EvidenceIndex[0][1].push_back(1428);
  fatSVCandidate.bp2EvidenceIndex[0][1].push_back(1480);
  fatSVCandidate.bp2EvidenceIndex[0][1].push_back(1507);

  std::vector<double> assemblyNoiseRate;
  assemblyNoiseRate.push_back(0.000002);  // 1st bam assembly noise rate
  assemblyNoiseRate.push_back(0.008);     // 2nd bam assembly noise rate

  // complex candidate is significant for 1st bam
  BOOST_REQUIRE(isAnyComplexCandidateSignalSignificant(1, fatSVCandidate, assemblyNoiseRate));

  assemblyNoiseRate.clear();
  assemblyNoiseRate.push_back(0.008);  // 1st bam assembly noise rate
  assemblyNoiseRate.push_back(0.008);  // 2nd bam assembly noise rate

  // complex candidate is not significant for any of the bams
  BOOST_REQUIRE(!isAnyComplexCandidateSignalSignificant(1, fatSVCandidate, assemblyNoiseRate));
}

// test the candidate's filtration state. This test verifies the following cases:
// 1. SEMI_MAPPED -  The candidate has ONLY evidence from read pairs where one read has MAPQ smaller than the
// cutoff
// (including MAPQ0)
// 2. SPANNING_LOW_SIGNAL - Candidates support spanning SV, but none of them is significant
// spanning candidate(Significant spanning SV candidate has been described in
// test_IsAnySpanningCandidateSignalSignificant)
// 3. COMPLEX_LOW_COUNT - When a complex SV doesn't satisfy minimum candidate count criteria
// 4. COMPLEX_LOW_SIGNAL - Candidates support complex SV, but none of them is significant
// complex candidate(Significant complex SV candidate has been described in
// test_IsAnyComplexCandidateSignalSignificant)
// 5. None - None of the above filtration state
BOOST_AUTO_TEST_CASE(test_IsFilterSingleJunctionCandidate)
{
  SVCandidate    svCandidate;
  FatSVCandidate fatSVCandidate1(svCandidate, 1u);

  fatSVCandidate1.bp1EvidenceIndex[0][0].push_back(3443);
  fatSVCandidate1.bp1EvidenceIndex[0][0].push_back(3452);
  fatSVCandidate1.bp1EvidenceIndex[0][0].push_back(3440);
  fatSVCandidate1.bp1EvidenceIndex[0][0].push_back(3489);

  fatSVCandidate1.bp2EvidenceIndex[0][0].push_back(1403);
  fatSVCandidate1.bp2EvidenceIndex[0][0].push_back(1428);
  fatSVCandidate1.bp2EvidenceIndex[0][0].push_back(1480);
  fatSVCandidate1.bp2EvidenceIndex[0][0].push_back(1507);

  std::vector<double> assemblyNoiseRate;
  assemblyNoiseRate.push_back(0.000002);
  std::vector<double> spanningNoiseRate;
  spanningNoiseRate.push_back(0.008);

  // test for SEMI_MAPPED candidate as filtration state
  BOOST_REQUIRE_EQUAL(
      isFilterSingleJunctionCandidate(false, spanningNoiseRate, assemblyNoiseRate, fatSVCandidate1, 1),
      SINGLE_JUNCTION_FILTER::SEMI_MAPPED);

  svCandidate.bp1.state = SVBreakendState::index_t::RIGHT_OPEN;
  svCandidate.bp2.state = SVBreakendState::index_t::LEFT_OPEN;
  svCandidate.bp1.lowresEvidence.add(0, 2);
  FatSVCandidate fatSVCandidate2(svCandidate, 1);
  fatSVCandidate2.bp1EvidenceIndex[0][0].push_back(3443);
  fatSVCandidate2.bp1EvidenceIndex[0][0].push_back(3452);
  fatSVCandidate2.bp1EvidenceIndex[0][0].push_back(3440);
  fatSVCandidate2.bp1EvidenceIndex[0][0].push_back(3489);

  fatSVCandidate2.bp2EvidenceIndex[0][0].push_back(1403);
  fatSVCandidate2.bp2EvidenceIndex[0][0].push_back(1428);
  fatSVCandidate2.bp2EvidenceIndex[0][0].push_back(1480);
  fatSVCandidate2.bp2EvidenceIndex[0][0].push_back(1507);

  // When none of the filtration state is satisfied
  spanningNoiseRate.clear();
  spanningNoiseRate.push_back(0.008);
  BOOST_REQUIRE_EQUAL(
      isFilterSingleJunctionCandidate(false, spanningNoiseRate, assemblyNoiseRate, fatSVCandidate2, 1),
      SINGLE_JUNCTION_FILTER::NONE);

  spanningNoiseRate.clear();
  spanningNoiseRate.push_back(0.1);
  // test for SPANNING_LOW_SIGNAL as filtration state
  BOOST_REQUIRE_EQUAL(
      isFilterSingleJunctionCandidate(false, spanningNoiseRate, assemblyNoiseRate, fatSVCandidate2, 1),
      SINGLE_JUNCTION_FILTER::SPANNING_LOW_SIGNAL);

  // test for COMPLEX_LOW_COUNT as filtration state
  svCandidate.bp1.state = SVBreakendState::index_t::COMPLEX;
  svCandidate.bp2.state = SVBreakendState::index_t::UNKNOWN;
  svCandidate.bp1.lowresEvidence.add(0, 2);
  FatSVCandidate fatSVCandidate3(svCandidate, 1);
  BOOST_REQUIRE_EQUAL(
      isFilterSingleJunctionCandidate(false, spanningNoiseRate, assemblyNoiseRate, fatSVCandidate3, 1),
      SINGLE_JUNCTION_FILTER::COMPLEX_LOW_COUNT);

  // test for COMPLEX_LOW_SIGNAL as filtration state
  svCandidate.bp1.state = SVBreakendState::index_t::COMPLEX;
  svCandidate.bp2.state = SVBreakendState::index_t::UNKNOWN;
  svCandidate.bp1.lowresEvidence.clear();
  svCandidate.bp1.lowresEvidence.add(2, 3);
  FatSVCandidate fatSVCandidate4(svCandidate, 1);
  BOOST_REQUIRE_EQUAL(
      isFilterSingleJunctionCandidate(false, spanningNoiseRate, assemblyNoiseRate, fatSVCandidate4, 1),
      SINGLE_JUNCTION_FILTER::COMPLEX_LOW_SIGNAL);
}

// test filters on all SV candidates. Following candidates will be filtered out
// 1. Semi Mapped
// 2. COMPLEX LOW COUNT
// 3. COMPLEX LOW SIGNAL
// This test also  checks the delaying process for filtering Spanning_Low_Signal candidates.
BOOST_AUTO_TEST_CASE(test_filterCandidates)
{
  SVCandidate    svCandidate1;
  FatSVCandidate fatSVCandidate1(svCandidate1, 1u);

  fatSVCandidate1.bp1EvidenceIndex[0][0].push_back(3443);
  fatSVCandidate1.bp1EvidenceIndex[0][0].push_back(3452);
  fatSVCandidate1.bp1EvidenceIndex[0][0].push_back(3440);
  fatSVCandidate1.bp1EvidenceIndex[0][0].push_back(3489);

  fatSVCandidate1.bp2EvidenceIndex[0][0].push_back(1403);
  fatSVCandidate1.bp2EvidenceIndex[0][0].push_back(1428);
  fatSVCandidate1.bp2EvidenceIndex[0][0].push_back(1480);
  fatSVCandidate1.bp2EvidenceIndex[0][0].push_back(1507);

  std::vector<double> assemblyNoiseRate;
  assemblyNoiseRate.push_back(0.000002);
  std::vector<double> spanningNoiseRate;
  spanningNoiseRate.push_back(0.008);

  std::vector<FatSVCandidate> svCandidates;
  // SEMI_MAPPED SV candidate. It should be filtered out.
  svCandidates.push_back(fatSVCandidate1);

  // spanning Low signal SV Cnadidate. It should not be filtered out.
  SVCandidate svCandidate2;
  svCandidate2.bp1.state = SVBreakendState::index_t::RIGHT_OPEN;
  svCandidate2.bp2.state = SVBreakendState::index_t::LEFT_OPEN;
  svCandidate2.bp1.lowresEvidence.add(0, 2);
  FatSVCandidate fatSVCandidate2(svCandidate2, 1);
  svCandidates.push_back(fatSVCandidate2);

  // COMPLEX LOW COUNT SV Candidate. It should be filtered out.
  SVCandidate svCandidate3;
  svCandidate3.bp1.state = SVBreakendState::index_t::COMPLEX;
  svCandidate3.bp2.state = SVBreakendState::index_t::UNKNOWN;
  svCandidate3.bp1.lowresEvidence.add(0, 2);
  FatSVCandidate fatSVCandidate3(svCandidate3, 1);
  svCandidates.push_back(fatSVCandidate3);

  // COMPLEX LOW SIGNAL. It should be filtered out.
  SVCandidate svCandidate4;
  svCandidate4.bp1.state = SVBreakendState::index_t::COMPLEX;
  svCandidate4.bp2.state = SVBreakendState::index_t::UNKNOWN;
  svCandidate4.bp1.lowresEvidence.clear();
  svCandidate4.bp1.lowresEvidence.add(2, 3);
  FatSVCandidate fatSVCandidate4(svCandidate4, 1);
  svCandidates.push_back(fatSVCandidate4);
  SVFinderStats stats;
  filterCandidates(false, spanningNoiseRate, assemblyNoiseRate, svCandidates, stats, 1);

  // check all the stats
  BOOST_REQUIRE_EQUAL(svCandidates.size(), 1);
  BOOST_REQUIRE_EQUAL(stats.ComplexLowCountFilter, 1);
  BOOST_REQUIRE_EQUAL(stats.ComplexLowSignalFilter, 1);
  BOOST_REQUIRE_EQUAL(stats.semiMappedFilter, 1);

  // check whether spanning Low signal sv candidate is there or not. It should not be filtered out.
  BOOST_REQUIRE_EQUAL(svCandidates[0].bp1.state, SVBreakendState::index_t::RIGHT_OPEN);
  BOOST_REQUIRE_EQUAL(svCandidates[0].bp2.state, SVBreakendState::index_t::LEFT_OPEN);
  // test whether spanning Low signal sv candidate is marked for a multi-junction evaluation.
  BOOST_REQUIRE(svCandidates[0].isSingleJunctionFilter);
}

// updateEvidenceIndex stores additional bam read index to decide if the candidate evidence
// is significant relative to background noise in the sample. The significance of SV candidate
// has been described in test_IsSpanningCandidateSignalSignificant and
// test_IsComplexCandidateSignalSignificant test cases. This unit test checks whether
// updateEvidenceIndex method stores read index correctly for different SV evidence. This test
// checks the read index based on nature of SV evidence provided by a single DNA/RNA fragment.
BOOST_AUTO_TEST_CASE(test_updateEvidenceIndex)
{
  bam_record bamRecord1;
  buildTestBamRecord(bamRecord1, 0, 200, 0, 210, 20, 15);
  bamRecord1.set_qname("Read-1");

  SVCandidateSetSequenceFragment fragment;
  SVObservation                  svObservation;
  // single source sv evidence
  svObservation.dnaFragmentSVEvidenceSource = SourceOfSVEvidenceInDNAFragment::READ1;
  fragment.read1.bamrec                     = bamRecord1;
  fragment.read1.readIndex                  = 1;  // setting the read index

  // check read index for Semi align evidence type
  svObservation.svEvidenceType = SVEvidenceType::SEMIALIGN;
  SVCandidate    svCandidate;
  FatSVCandidate fatSVCandidate(svCandidate, 1u);
  updateEvidenceIndex(fragment, svObservation, fatSVCandidate, 0);
  BOOST_REQUIRE_EQUAL(fatSVCandidate.bp1EvidenceIndex[SVEvidenceType::SEMIALIGN][0][0], 1);

  // check read index for SPLIT_ALIGN align evidence type
  svObservation.svEvidenceType = SVEvidenceType::SPLIT_ALIGN;
  updateEvidenceIndex(fragment, svObservation, fatSVCandidate, 0);
  BOOST_REQUIRE_EQUAL(fatSVCandidate.bp1EvidenceIndex[SVEvidenceType::SPLIT_ALIGN][0].size(), 1);
  BOOST_REQUIRE_EQUAL(fatSVCandidate.bp2EvidenceIndex[SVEvidenceType::SPLIT_ALIGN][0].size(), 0);
  BOOST_REQUIRE_EQUAL(fatSVCandidate.bp1EvidenceIndex[SVEvidenceType::SPLIT_ALIGN][0][0], 1);

  // Adding supplementary read
  bam_record supplementSASplitRead;
  buildTestBamRecord(supplementSASplitRead);
  addSupplementaryAlignmentEvidence(supplementSASplitRead);

  SVCandidateSetRead svCandidateSetRead;
  svCandidateSetRead.bamrec = supplementSASplitRead;
  fragment.read1Supplemental.push_back(svCandidateSetRead);
  updateEvidenceIndex(fragment, svObservation, fatSVCandidate, 0);
  BOOST_REQUIRE_EQUAL(fatSVCandidate.bp1EvidenceIndex[SVEvidenceType::SPLIT_ALIGN][0].size(), 2);
  BOOST_REQUIRE_EQUAL(fatSVCandidate.bp2EvidenceIndex[SVEvidenceType::SPLIT_ALIGN][0].size(), 1);
  BOOST_REQUIRE_EQUAL(fatSVCandidate.bp1EvidenceIndex[SVEvidenceType::SPLIT_ALIGN][0][1], 1);
  BOOST_REQUIRE_EQUAL(fatSVCandidate.bp2EvidenceIndex[SVEvidenceType::SPLIT_ALIGN][0][0], 0);

  // check read index for PAIR evidence type
  bam_record bamRecord2;
  buildTestBamRecord(bamRecord2, 0, 210, 0, 220, 20, 15);
  bamRecord2.set_qname("Read-2");

  // SV evidence source is Pair Reads
  svObservation.dnaFragmentSVEvidenceSource = SourceOfSVEvidenceInDNAFragment::READ_PAIR;
  fragment.read1Supplemental.clear();
  fragment.read2.bamrec        = bamRecord2;
  fragment.read2.readIndex     = 2;
  svObservation.svEvidenceType = SVEvidenceType::PAIR;
  updateEvidenceIndex(fragment, svObservation, fatSVCandidate, 0);
  BOOST_REQUIRE_EQUAL(fatSVCandidate.bp1EvidenceIndex[SVEvidenceType::PAIR][0].size(), 1);
  BOOST_REQUIRE_EQUAL(fatSVCandidate.bp2EvidenceIndex[SVEvidenceType::PAIR][0].size(), 1);
  BOOST_REQUIRE_EQUAL(fatSVCandidate.bp1EvidenceIndex[SVEvidenceType::PAIR][0][0], 1);
  BOOST_REQUIRE_EQUAL(fatSVCandidate.bp2EvidenceIndex[SVEvidenceType::PAIR][0][0], 2);
}

/// Check whether any SVs can intersect each other
/// Two SVCandidates intersect if their breakend regions overlap in the same direction.
/// The test case examples are illustrated below, where the sequentially
/// intersecting candidates are (1, 2, 5, 6).
///
/// fatSVCandidate Index and Schematic:
/// 1: >>>bp1>>>>-------------------------<<bp2<<<<< [bp1(10,100) & bp2(1000,1100)]
/// 2:   >>>>bp1>>>--------------------------<<bp2<<<<< [bp1(50,120) & bp2(1050,1150)]
/// 3:               >>>bp1>>>-------------------------------<<<bp2<<<<< [bp1(200,300) & bp2(1400,1500)]
/// 4:   <<<bp1<<<<---------------------<<<<<<<bp2<<<<< [bp1(50,120) & bp2(990,1050)]
/// 5:        >>>bp1>>>>--------------------------<<bp2<<<<< [bp1(110,150) & bp2(1110, 1250)]
/// 6:               >>>bp1>>>>-------------------------<<bp2<<<<< [bp1(140, 200) & bp2(1240, 1300)]
BOOST_AUTO_TEST_CASE(test_consolidateOverlap)
{
  SVCandidate svCandidate1;
  svCandidate1.bp1.interval = GenomeInterval(0, 10, 100);
  svCandidate1.bp2.interval = GenomeInterval(0, 1000, 1100);
  svCandidate1.bp1.state    = SVBreakendState::RIGHT_OPEN;
  svCandidate1.bp2.state    = SVBreakendState::LEFT_OPEN;

  SVCandidate svCandidate2;
  svCandidate2.bp1.interval = GenomeInterval(0, 50, 120);
  svCandidate2.bp2.interval = GenomeInterval(0, 1050, 1150);
  svCandidate2.bp1.state    = SVBreakendState::RIGHT_OPEN;
  svCandidate2.bp2.state    = SVBreakendState::LEFT_OPEN;

  SVCandidate svCandidate3;
  svCandidate3.bp1.interval = GenomeInterval(0, 200, 300);
  svCandidate3.bp2.interval = GenomeInterval(0, 1400, 1500);
  svCandidate3.bp1.state    = SVBreakendState::RIGHT_OPEN;
  svCandidate3.bp2.state    = SVBreakendState::LEFT_OPEN;

  SVCandidate svCandidate4;
  svCandidate3.bp1.interval = GenomeInterval(0, 50, 120);
  svCandidate3.bp2.interval = GenomeInterval(0, 990, 1050);
  svCandidate3.bp1.state    = SVBreakendState::LEFT_OPEN;
  svCandidate3.bp2.state    = SVBreakendState::LEFT_OPEN;

  FatSVCandidate fatSVCandidate1(svCandidate1, 1u);
  FatSVCandidate fatSVCandidate2(svCandidate2, 1u);
  FatSVCandidate fatSVCandidate3(svCandidate3, 1u);
  FatSVCandidate fatSVCandidate4(svCandidate4, 1u);

  std::vector<FatSVCandidate> fatSVs;
  fatSVs.push_back(fatSVCandidate1);
  fatSVs.push_back(fatSVCandidate3);
  fatSVs.push_back(fatSVCandidate4);

  const bam_header_info bamHeader(buildTestBamHeader());
  std::string           queryseq1 = "AGCTGACTGATCGATTTTTTACGTAGAGGAGCTTTGACGTATGAGCCTGATATGAGCCTG";
  std::string           queryseq2 = "TGACGTATGAGCCTGATATGAGCCT";
  bam_record            bamRecord1;
  buildTestBamRecord(bamRecord1, 0, 200, 0, 300, 60, 15, "60M", queryseq1);
  bamRecord1.set_qname("Read-1");
  bamRecord1.toggle_is_first();

  bam_record bamRecord2;
  buildTestBamRecord(bamRecord2, 0, 300, 0, 200, 25, 15, "25M", queryseq2);
  bamRecord2.toggle_is_second();
  bamRecord2.toggle_is_fwd_strand();
  bamRecord2.toggle_is_mate_fwd_strand();
  bamRecord2.set_qname("Read-1");
  SVCandidateSetData                         candidateSetData;
  SVCandidateSetSequenceFragmentSampleGroup& group = candidateSetData.getDataGroup(0);
  group.add(bamHeader, bamRecord1, false, true, true);
  group.add(bamHeader, bamRecord2, false, true, true);
  SVSequenceFragmentAssociation association(0, SVEvidenceType::PAIR);
  group.begin().                operator*().svLink.push_back(association);
  // candidate-1, candidate-3 and candidate-4 are not overlapping
  // So consolidated size = 3.
  consolidateOverlap(1, candidateSetData, fatSVs);
  BOOST_REQUIRE_EQUAL(fatSVs.size(), 3);

  fatSVs.clear();
  SVCandidate svCandidate5;
  svCandidate5.bp1.interval = GenomeInterval(0, 110, 150);
  svCandidate5.bp2.interval = GenomeInterval(0, 1110, 1250);
  svCandidate5.bp1.state    = SVBreakendState::RIGHT_OPEN;
  svCandidate5.bp2.state    = SVBreakendState::LEFT_OPEN;
  FatSVCandidate fatSVCandidate5(svCandidate5, 1u);
  SVCandidate    svCandidate6;
  svCandidate6.bp1.interval = GenomeInterval(0, 140, 200);
  svCandidate6.bp2.interval = GenomeInterval(0, 1240, 1300);
  svCandidate6.bp1.state    = SVBreakendState::RIGHT_OPEN;
  svCandidate6.bp2.state    = SVBreakendState::LEFT_OPEN;
  FatSVCandidate fatSVCandidate6(svCandidate6, 1u);
  fatSVs.push_back(fatSVCandidate1);
  fatSVs.push_back(fatSVCandidate2);
  fatSVs.push_back(fatSVCandidate4);
  fatSVs.push_back(fatSVCandidate5);
  fatSVs.push_back(fatSVCandidate6);
  // Candidate-1 overlaps with candidate2. Resultant candidate overlaps with
  // Candidate-5, agian resultant candidate overlaps with candidate-6
  // So consolidated size = 2.
  consolidateOverlap(1, candidateSetData, fatSVs);
  BOOST_REQUIRE_EQUAL(fatSVs.size(), 2);
}

// Test the following cases:
// 1. When number of edges in a SV locus less than min edge count, no sv will be generated for this locus.
// 2. Edge should be bidirectional that means if there is an edge from node-1 to node-2, there should be an
//    edge from node-2 to node-1
// 3. If min edge criteria is satisfied, then for a bam read it will return a sv candidate.
BOOST_AUTO_TEST_CASE(test_SVCandidates)
{
  const bam_header_info           bamHeader(buildTestBamHeader());
  std::unique_ptr<SVLocusScanner> scanner(buildTestSVLocusScanner(bamHeader));
  const std::string               referenceFilename = getTestReferenceFilename();

  // 1. Local Evidence = GenomeInterval: 0:[80,120) (Read posotion is 100 and 20 bp padding on each side)
  // 2. Remote Evidence = GenomeInterval: 0:[279,319) (According to SA tag other segment position is 299 and
  // 20bp
  bam_record supplementSASplitRead;
  buildTestBamRecord(supplementSASplitRead);                 // pos=100, matePos=200 and fragmentSize=100
  addSupplementaryAlignmentEvidence(supplementSASplitRead);  // SA tag : chrFoo,300,-,54H22M,50,0;
  std::vector<bam_record> readsToAdd;
  readsToAdd.push_back(supplementSASplitRead);
  const std::shared_ptr<BamFilenameMaker> bamFileNameMaker(new BamFilenameMaker());
  const std::string&                      bamFileName(bamFileNameMaker.get()->getFilename());
  buildTestBamFile(bamHeader, readsToAdd, bamFileName);
  // Build chromosome depth file
  const std::shared_ptr<TestChromosomeDepthFileMaker> depthFileMaker(new TestChromosomeDepthFileMaker());
  const std::string                                   depthFileName(depthFileMaker.operator*().getFilename());
  buildTestChromosomeDepthFile(depthFileName);

  // Create options
  TestFilenameMaker testFilenameMaker;
  TestFilenameMaker fileNameMaker0;
  TestFilenameMaker fileNameMaker1;
  GSCOptions        options;
  options.edgeRuntimeFilename             = fileNameMaker0.getFilename();
  options.edgeStatsFilename               = fileNameMaker1.getFilename();
  options.referenceFilename               = referenceFilename;
  options.alignFileOpt.alignmentFilenames = {bamFileName};
  options.alignFileOpt.isAlignmentTumor   = {false};
  options.graphFilename                   = testFilenameMaker.getFilename();
  ;
  options.chromDepthFilename = depthFileName;

  // Various other data structures shared by all tests:
  SVLocus locus1;
  locusAddPair(locus1, 0, 80, 120, 0, 279, 319);
  SVLocus locus2;
  locusAddPair(locus2, 0, 279, 319, 0, 80, 120);
  SVLocusSetOptions sopt;

  static const bool isSkipLocusSetIndexCreation(true);

  auto                edgeTrackerPtr(std::make_shared<EdgeRuntimeTracker>(options.edgeRuntimeFilename));
  GSCEdgeStatsManager edgeStatMan;

  SVCandidateSetData       svData;
  std::vector<SVCandidate> svs;
  EdgeInfo                 edgeInfo;
  edgeInfo.nodeIndex1 = 0;
  edgeInfo.nodeIndex2 = 1;

  // Create 1st SVLocus set
  {
    sopt.minMergeEdgeObservations = 10;
    SVLocusSet set(sopt, bamHeader, {bamFileName});
    set.merge(locus1);
    set.merge(locus2);
    set.checkState(true, true);
    // Serialize SV locus graph
    set.save(options.graphFilename.c_str());

    // Deserialize SV locus graph
    const SVLocusSet cset(options.graphFilename.c_str(), isSkipLocusSetIndexCreation);
    SVFinder         finder(
                options,
                scanner.operator*(),
                cset.getBamHeader(),
                cset.getAllSampleReadCounts(),
                edgeTrackerPtr,
                edgeStatMan);
    finder.findCandidateSV(cset, edgeInfo, svData, svs);
    // Min edge criteria is not satisfied as minimum
    // number of edges required is 10
    BOOST_REQUIRE_EQUAL(svData.getDataGroup(0).size(), 0);
    BOOST_REQUIRE_EQUAL(svs.size(), 0);
  }

  // Designed the Case-2 where there is an edge from locus node [80,120) to
  // locus node [279,319) but there is no edge from locus node [279,319) to
  // locus node [80,120). So number of SV candidates are 0. Here minimum number
  // of edges required is 1.
  {
    sopt.minMergeEdgeObservations = 1;
    SVLocusSet set(sopt, bamHeader, {bamFileName});
    // One edge is created.
    set.merge(locus1);
    set.checkState(true, true);
    // Serialize SV locus graph
    set.save(options.graphFilename.c_str());

    // Deserialize SV locus graph
    const SVLocusSet cset(options.graphFilename.c_str(), isSkipLocusSetIndexCreation);
    SVFinder         finder2(
                options,
                scanner.operator*(),
                cset.getBamHeader(),
                cset.getAllSampleReadCounts(),
                edgeTrackerPtr,
                edgeStatMan);
    finder2.findCandidateSV(cset, edgeInfo, svData, svs);
    BOOST_REQUIRE_EQUAL(svData.getDataGroup(0).size(), 0);
    BOOST_REQUIRE_EQUAL(svs.size(), 0);
  }

  // Designed the Case-3 where there is an edge from locus node [80,120) to
  // locus node [279,319) and also there is an edge from locus node [279,319) to
  // locus node [80,120). Here minimum number of edges required is 1.
  {
    sopt.minMergeEdgeObservations = 1;
    SVLocusSet set(sopt, bamHeader, {bamFileName});
    set.merge(locus1);
    set.merge(locus2);
    set.checkState(true, true);
    // Serialize SV locus graph
    set.save(options.graphFilename.c_str());

    // Deserialize SV locus graph
    const SVLocusSet cset(options.graphFilename.c_str(), isSkipLocusSetIndexCreation);
    SVFinder         finder(
                options,
                scanner.operator*(),
                cset.getBamHeader(),
                cset.getAllSampleReadCounts(),
                edgeTrackerPtr,
                edgeStatMan);
    finder.findCandidateSV(cset, edgeInfo, svData, svs);
    BOOST_REQUIRE_EQUAL(svData.getDataGroup(0).size(), 1);
    BOOST_REQUIRE_EQUAL(svs.size(), 1);
  }
}

BOOST_AUTO_TEST_SUITE_END()
