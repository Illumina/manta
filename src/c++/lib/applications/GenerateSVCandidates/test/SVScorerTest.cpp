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
#include "test/testAlignmentDataUtil.hpp"
#include "test/testFileMakers.hpp"
#include "test/testSVLocusScanner.hpp"

#include "SVScorer.cpp"
#include "SVScorer.hpp"

/// TestSVScorer is a friend of SVScorer. So that can access private
/// method of SVScorer
struct TestSVScorer {
  void getBreakendMaxMappedDepthAndMQ0(
      SVScorer&         scorer,
      const bool        isTumorOnly,
      const bool        isMaxDepth,
      const double      cutoffDepth,
      const SVBreakend& breakend,
      unsigned&         maxDepth,
      float&            MQ0Frac)
  {
    scorer.getBreakendMaxMappedDepthAndMQ0(isTumorOnly, isMaxDepth, cutoffDepth, breakend, maxDepth, MQ0Frac);
  }

  void computeAllScoreModels(
      SVScorer&                            scorer,
      const bool                           isSomatic,
      const bool                           isTumorOnly,
      const std::vector<JunctionCallInfo>& junctionData,
      SVModelScoreInfo&                    modelScoreInfo)
  {
    if (isSomatic) {
      scorer._sampleCount        = 2;
      scorer._diploidSampleCount = 1;
    }
    scorer.computeAllScoreModels(isSomatic, isTumorOnly, junctionData, modelScoreInfo);
  }

  void setSampleCount(SVScorer& scorer, unsigned sampleCount, unsigned diploidSampleCount)
  {
    scorer._sampleCount        = sampleCount;
    scorer._diploidSampleCount = diploidSampleCount;
  }
};

BOOST_AUTO_TEST_SUITE(SVScorer_test_suite)

// Create Temporary bam streams of a bam file which contains
// twenty one bam records. Ou of twenty one bam records, 18 are
// with mapping quality 0.
struct BamStream {
  BamStream()
  {
    const bam_header_info bamHeader(buildTestBamHeader());

    std::string querySeq = "TCTATCACCCATTTTACCACTCACGGGAGCTCTCC";
    bam_record  bamRecord1;
    buildTestBamRecord(bamRecord1, 0, 9, 0, 100, 35, 15, "35M", querySeq, 200);
    bamRecord1.set_qname("bamRecord1");
    // small Fragment length read
    bam_record bamRecord2;
    buildTestBamRecord(bamRecord2, 0, 109, 0, 125, 35, 15, "35M", querySeq, 49);
    bamRecord2.set_qname("bamRecord2");
    bam_record bamRecord3;
    buildTestBamRecord(bamRecord3, 0, 109, 0, 180, 35, 15, "35M", querySeq, 100);
    bamRecord3.set_qname("bamRecord3");
    readsToAdd.push_back(bamRecord1);
    readsToAdd.push_back(bamRecord2);
    readsToAdd.push_back(bamRecord3);
    for (unsigned i(0); i < 18; i++) {
      bam_record bamRecord;
      buildTestBamRecord(bamRecord, 0, 109, 0, 180, 35, 0, "35M", querySeq, 100);
      std::string name = "bamRecord" + std::to_string(i + 4);
      bamRecord.set_qname(name.c_str());
      readsToAdd.push_back(bamRecord);
    }
    bamFilename = _bamFilename();
    buildTestBamFile(bamHeader, readsToAdd, bamFilename);

    const std::string                          referenceFilename = getTestReferenceFilename();
    std::vector<std::string>                   bamFilenames      = {bamFilename};
    std::vector<std::shared_ptr<bam_streamer>> bamStreams;
    openBamStreams(referenceFilename, bamFilenames, bamStreams);
    bamStream = bamStreams[0];
  }

  std::shared_ptr<bam_streamer> bamStream;
  std::vector<bam_record>       readsToAdd;
  std::string                   bamFilename;

private:
  const std::string&     _bamFilename() const { return _bamFilenameMaker.getFilename(); }
  const BamFilenameMaker _bamFilenameMaker;
};

BOOST_FIXTURE_TEST_SUITE(SVScorer_test_suite, BamStream)

// test depth on each location i.e. number of
// read bases overlap in a location.
// We have taken here small fragment legth read as
// here we are testing the utility.
BOOST_AUTO_TEST_CASE(test_addReadToDepthEst)
{
  // read length = 15
  bam_record bamRecord1;
  buildTestBamRecord(bamRecord1, 0, 200, 0, 210, 15, 15, "15M");
  bamRecord1.set_qname("Read-1");

  // Read length = 15
  bam_record bamRecord2;
  buildTestBamRecord(bamRecord2, 0, 210, 0, 220, 15, 15, "15M");
  bamRecord2.set_qname("Read-2");

  std::vector<unsigned> depth(30);
  addReadToDepthEst(bamRecord1, 200, depth);
  addReadToDepthEst(bamRecord2, 200, depth);

  // test the coverage
  for (unsigned i = 0; i < 25; i++) {
    // From 210 to 214 i.e. total 5 locations are common in both the reads
    if (i >= 10 && i <= 14)
      BOOST_REQUIRE_EQUAL(depth[i], 2u);
    else
      BOOST_REQUIRE_EQUAL(depth[i], 1u);
  }
}

// Test the following simple utility mathematics
// lnToProb(a, b): a = exp(a-b)/(1 + exp(a-b))
//                 b = 1/(1 + exp(a-b))
BOOST_AUTO_TEST_CASE(test_lnToProb)
{
  float               lower(5);
  float               higher(5);
  static const double eps(0.00000001);
  // exp(a-b) = 1. so a = 1/2 and b = 1/2
  lnToProb(lower, higher);
  BOOST_REQUIRE_CLOSE(lower, 0.5, eps);
  BOOST_REQUIRE_CLOSE(higher, 0.5, eps);
}

// Test the following cases
// 1. If forced support is not allowed, a read should support split read evidence either at
//    breakpoint-1 or breakpoint-2 and calculate lnlhood according to case-3 and case-4.
// 2. If forced support is allowed, then even if the read is not supporting the breakpoints,
//    altlnlhood (maximum of lnlhood of BP1 and BP2)and refLnLhood (maximum of lnlhood of
//    BP1 and BP2) are still calculated.
// 3. If forced support is not allowed Value of altlnlhood should be calculated as follows:
//        altlnlhood = maximum of bp1 lnlhood and bp2 lnlhood of alt allele if
//        both the breakpoints support split evidence, otherwise corresponding breakpoint
//        lnlhood will be taken.
// 4. If forced support is not allowed, refLnLhood should be maximum of bp1 lnlhood and bp2
//    lnlhood of Ref allele if both the breakpoints support split evidence, otherwise corresponding breakpoint
//    lnlhood will be taken.
BOOST_AUTO_TEST_CASE(test_getSampleSplitReadLnLhood)
{
  float               refSplitLnLhoodBP1(-0.5);
  float               refSplitLnLhoodBP2(-0.8);
  float               altSplitLnLhoodBP1(-0.2);
  float               altSplitLnLhoodBP2(-0.4);
  static const double eps(0.00000001);
  float               refLnLhood;
  float               altLnLhood;
  SVFragmentEvidence  fragmentEvidence;

  // Both the breakpoints are not supporting split read evidence
  // for ref and alt alllele, but forced support is allowed.
  fragmentEvidence.alt.bp1.read1.isSplitSupport = false;
  fragmentEvidence.alt.bp1.read1.splitLnLhood   = altSplitLnLhoodBP1;
  fragmentEvidence.alt.bp2.read1.isSplitSupport = false;
  fragmentEvidence.alt.bp2.read1.splitLnLhood   = altSplitLnLhoodBP2;
  fragmentEvidence.ref.bp1.read1.isSplitSupport = false;
  fragmentEvidence.ref.bp1.read1.splitLnLhood   = refSplitLnLhoodBP1;
  fragmentEvidence.ref.bp2.read1.isSplitSupport = false;
  fragmentEvidence.ref.bp2.read1.splitLnLhood   = refSplitLnLhoodBP2;
  // refLnLhood = maximum of (refSplitLnLhoodBP1, refSplitLnLhoodBP2)
  // altLnLhood = maximum of (altSplitLnLhoodBP1, altSplitLnLhoodBP2)
  BOOST_REQUIRE(getSampleSplitReadLnLhood(fragmentEvidence, true, refLnLhood, altLnLhood, true));
  BOOST_REQUIRE_CLOSE(refLnLhood, refSplitLnLhoodBP1, eps);
  BOOST_REQUIRE_CLOSE(altLnLhood, altSplitLnLhoodBP1, eps);

  // BP1 and BP2 of alt allele are supporting split read evidence.
  fragmentEvidence.alt.bp1.read1.isSplitSupport = true;
  fragmentEvidence.alt.bp1.read1.splitLnLhood   = altSplitLnLhoodBP1;
  fragmentEvidence.alt.bp2.read1.isSplitSupport = true;
  fragmentEvidence.alt.bp2.read1.splitLnLhood   = altSplitLnLhoodBP2;
  fragmentEvidence.ref.bp1.read1.splitLnLhood   = refSplitLnLhoodBP1;
  fragmentEvidence.ref.bp2.read1.splitLnLhood   = refSplitLnLhoodBP2;
  BOOST_REQUIRE(getSampleSplitReadLnLhood(fragmentEvidence, true, refLnLhood, altLnLhood));
  // maximum of (refSplitLnLhoodBP1, refSplitLnLhoodBP2)
  BOOST_REQUIRE_CLOSE(refLnLhood, refSplitLnLhoodBP1, eps);
  // maximum of (altSplitLnLhoodBP1, altSplitLnLhoodBP2)
  BOOST_REQUIRE_CLOSE(altLnLhood, altSplitLnLhoodBP1, eps);

  // Read-1 is supporting split evidence of alt allele only on BP1
  fragmentEvidence.alt.bp1.read1.isSplitSupport = true;
  fragmentEvidence.ref.bp2.read1.isSplitSupport = false;
  BOOST_REQUIRE(getSampleSplitReadLnLhood(fragmentEvidence, true, refLnLhood, altLnLhood));
  // As split support evidence is present only on BP1, then lnlhood of BP1 will be taken.
  BOOST_REQUIRE_CLOSE(refLnLhood, refSplitLnLhoodBP1, eps);
  BOOST_REQUIRE_CLOSE(altLnLhood, altSplitLnLhoodBP1, eps);

  // Read-1 is supporting split evidence of alt allele only on BP2
  fragmentEvidence.alt.bp1.read1.isSplitSupport = false;
  fragmentEvidence.alt.bp2.read1.isSplitSupport = true;
  BOOST_REQUIRE(getSampleSplitReadLnLhood(fragmentEvidence, true, refLnLhood, altLnLhood));
  // As split support evidence is present only on BP2, then lnlhood of BP2 will be taken.
  BOOST_REQUIRE_CLOSE(refLnLhood, refSplitLnLhoodBP2, eps);
  BOOST_REQUIRE_CLOSE(altLnLhood, altSplitLnLhoodBP2, eps);

  // Both the breakpoints are not supporting split read evidence
  // for ref and alt alllele. API returns the default value of refLnLhood
  // and altLnLhood which is 1.
  fragmentEvidence.alt.bp1.read1.isSplitSupport = false;
  fragmentEvidence.alt.bp1.read1.splitLnLhood   = altSplitLnLhoodBP1;
  fragmentEvidence.alt.bp2.read1.isSplitSupport = false;
  fragmentEvidence.alt.bp2.read1.splitLnLhood   = altSplitLnLhoodBP2;
  fragmentEvidence.ref.bp1.read1.isSplitSupport = false;
  fragmentEvidence.ref.bp1.read1.splitLnLhood   = refSplitLnLhoodBP1;
  fragmentEvidence.ref.bp2.read1.isSplitSupport = false;
  fragmentEvidence.ref.bp2.read1.splitLnLhood   = refSplitLnLhoodBP2;
  BOOST_REQUIRE(!getSampleSplitReadLnLhood(fragmentEvidence, true, refLnLhood, altLnLhood));
  BOOST_REQUIRE_EQUAL(refLnLhood, 1.);
  BOOST_REQUIRE_EQUAL(altLnLhood, 1.);
}

// Test the following cases:
// 1. If altLnLhood is greater than refLnLhood and normalized probability
//    (explained in test_lnToProb) of altLnLhood is greater than splitSupportProb(0.999f),
//    then it is a confident split read count for alt allele.
// 2. If altLnLhood is greater than refLnLhood but normalized probability
//    (explained in test_lnToProb) of altLnLhood is not greater than splitSupportProb(0.999f),
//    then it is not a confident split read count for alt allele.
// 3. If refLnLhood is greater than altLnLhood and normalized probability
//    of altLnLhood is greater than splitSupportProb(0.999f), then it is
//    a confident split read count for ref allele.
// 4. For ref allele, if above case-3 is satisfied, then track the confident
//    split read count of BP1 and BP2.
// Normalized probability is calculated as,
// lnToProb(a, b): a = exp(a-b)/(1 + exp(a-b))
//                 b = 1/(1 + exp(a-b))
// where a and b are lnlhood values.
BOOST_AUTO_TEST_CASE(test_addConservativeSplitReadSupport)
{
  // Designed the case-1
  SVSampleInfo       sampleInfo1;
  SVFragmentEvidence fragmentEvidence;
  fragmentEvidence.alt.bp1.read1.isSplitSupport = true;
  fragmentEvidence.alt.bp1.read1.splitLnLhood   = -7.9;
  fragmentEvidence.alt.bp2.read1.splitLnLhood   = -8.9;
  fragmentEvidence.ref.bp1.read1.splitLnLhood   = -17.2;
  fragmentEvidence.ref.bp2.read1.splitLnLhood   = -18.9;
  // So altLnLhood = max(bp1 splitLnLhood of alt, bp2 splitLnLhood of alt) = -7.9
  // So refLnLhood = max(bp1 splitLnLhood of ref, bp2 splitLnLhood of ref) = -17.2
  // Here, altLnLhood > refLnLhood. So normalized probability(based on test_lnToProb)
  // of altLnLhood  = 1/exp(refLnLhood - altLnLhood)
  //                = 0.999909 which is greater than 0.999f
  // So it is a confident split read evidence for alt allele.
  addConservativeSplitReadSupport(fragmentEvidence, true, sampleInfo1);
  BOOST_REQUIRE_EQUAL(sampleInfo1.alt.confidentSplitReadCount, 1);

  // Designed the case-2
  fragmentEvidence.alt.bp1.read1.splitLnLhood = -7.9;
  fragmentEvidence.alt.bp2.read1.splitLnLhood = -8.9;
  fragmentEvidence.ref.bp1.read1.splitLnLhood = -10.2;
  fragmentEvidence.ref.bp2.read1.splitLnLhood = -18.9;
  SVSampleInfo sampleInfo2;
  // So altLnLhood = max(bp1 splitLnLhood of alt, bp2 splitLnLhood of alt) = -7.9
  // So refLnLhood = max(bp1 splitLnLhood of ref, bp2 splitLnLhood of ref) = -10.2
  // Here, altLnLhood > refLnLhood. So normalized probability(based on test_lnToProb)
  // of altLnLhood  = 1/exp(refLnLhood - altLnLhood)
  //                = 0.908877 which is less than 0.999f
  // So it is not a confident split read evidence for alt allele.
  addConservativeSplitReadSupport(fragmentEvidence, true, sampleInfo2);
  BOOST_REQUIRE_EQUAL(sampleInfo2.alt.confidentSplitReadCount, 0);

  // Designed the case-3
  SVSampleInfo sampleInfo3;
  fragmentEvidence.ref.bp1.read1.splitLnLhood   = -7.9;
  fragmentEvidence.ref.bp2.read1.splitLnLhood   = -8.9;
  fragmentEvidence.alt.bp1.read1.splitLnLhood   = -17.2;
  fragmentEvidence.alt.bp2.read1.splitLnLhood   = -18.9;
  fragmentEvidence.ref.bp1.read1.isSplitSupport = true;
  // So altLnLhood = max(bp1 splitLnLhood of alt, bp2 splitLnLhood of alt) = -17.2
  // So refLnLhood = max(bp1 splitLnLhood of ref, bp2 splitLnLhood of ref) = -7.9
  // Here, refLnLhood > altLnLhood. So normalized probability(based on test_lnToProb)
  // of refLnLhood  = 1/exp(altLnLhood - refLnLhood)
  //                = 0.999909 which is greater than 0.999f
  // So it is a confident split read evidence for ref allele.
  addConservativeSplitReadSupport(fragmentEvidence, true, sampleInfo3);
  BOOST_REQUIRE_EQUAL(sampleInfo3.ref.confidentSplitReadCount, 1);
  // Designed the case-4 where ref allele supports split evidence on BP1.
  BOOST_REQUIRE_EQUAL(sampleInfo3.ref.confidentSplitReadAndPairCountRefBp1, 1);

  // Designed the case-3
  SVSampleInfo sampleInfo4;
  fragmentEvidence.ref.bp1.read1.splitLnLhood   = -8.9;
  fragmentEvidence.ref.bp2.read1.splitLnLhood   = -7.9;
  fragmentEvidence.alt.bp1.read1.splitLnLhood   = -17.2;
  fragmentEvidence.alt.bp2.read1.splitLnLhood   = -18.9;
  fragmentEvidence.ref.bp1.read1.isSplitSupport = false;
  fragmentEvidence.ref.bp2.read1.isSplitSupport = true;
  // So altLnLhood = max(bp1 splitLnLhood of alt, bp2 splitLnLhood of alt) = -17.2
  // So refLnLhood = max(bp1 splitLnLhood of ref, bp2 splitLnLhood of ref) = -7.9
  // Here, refLnLhood > altLnLhood. So normalized probability(based on test_lnToProb)
  // of refLnLhood  = 1/exp(altLnLhood - refLnLhood)
  //                = 0.999909 which is greater than 0.999f
  // So it is a confident split read evidence for ref allele.
  addConservativeSplitReadSupport(fragmentEvidence, true, sampleInfo4);
  BOOST_REQUIRE_EQUAL(sampleInfo4.ref.confidentSplitReadCount, 1);
  // Designed the case-4 where ref allele supports split evidence on BP2.
  BOOST_REQUIRE_EQUAL(sampleInfo4.ref.confidentSplitReadAndPairCountRefBp2, 1);
}

// Test the following cases:
// 1. If a fragment supports only breakpoint-1, return fragment size
//    probability of breakpoint-1
// 2. If a fragment supports only breakpoint-2, return fragment size
//    probability of breakpoint-2
// 3. If a fragment supports both breakpoint-1 and breakpoint-2, return maximum
//    of fragment size probability of breakpoint-1 and breakpoint-2
// 4. If a fragment does not support breakpoint-1 or breakpoint-2, return 0.
BOOST_AUTO_TEST_CASE(test_getSpanningPairAlleleLhood)
{
  float               fragLengthProb1(0.4);
  float               fragLengthProb2(0.6);
  static const double eps(0.00000001);

  // Designed the case-1 where a fragment supports allele only in BP1
  SVFragmentEvidenceAllele fragmentEvidenceAllele;
  fragmentEvidenceAllele.bp1.isFragmentSupport = true;
  fragmentEvidenceAllele.bp1.fragLengthProb    = fragLengthProb1;
  BOOST_REQUIRE_CLOSE(getSpanningPairAlleleLhood(fragmentEvidenceAllele), fragLengthProb1, eps);

  // Designed the case-2 where a fragment supports allele only in BP2
  fragmentEvidenceAllele.bp1.isFragmentSupport = false;
  fragmentEvidenceAllele.bp2.isFragmentSupport = true;
  fragmentEvidenceAllele.bp2.fragLengthProb    = fragLengthProb2;
  BOOST_REQUIRE_CLOSE(getSpanningPairAlleleLhood(fragmentEvidenceAllele), fragLengthProb2, eps);

  // Designed the case-3 where a fragment supports allele in both BP1 and BP2
  fragmentEvidenceAllele.bp1.isFragmentSupport = true;
  BOOST_REQUIRE_CLOSE(getSpanningPairAlleleLhood(fragmentEvidenceAllele), fragLengthProb2, eps);

  // Designed the case-4 where a fragment doesn't support allele either of the breakpoints.
  fragmentEvidenceAllele.bp1.isFragmentSupport = false;
  fragmentEvidenceAllele.bp2.isFragmentSupport = false;
  BOOST_REQUIRE_EQUAL(getSpanningPairAlleleLhood(fragmentEvidenceAllele), 0);
}

// Test Spanning pair count based on whether a fragment supports
// breakpoint-1(BP1) or breakpoint-2(BP2) for alt and ref allele.
BOOST_AUTO_TEST_CASE(test_addSpanningPairSupport)
{
  // Fragment does not support BP1 or BP2 for both alt and ref allele
  SVFragmentEvidence fragmentEvidence;
  SVSampleInfo       sampleInfo1;
  addSpanningPairSupport(fragmentEvidence, sampleInfo1);
  BOOST_REQUIRE_EQUAL(sampleInfo1.alt.spanningPairCount, 0);
  BOOST_REQUIRE_EQUAL(sampleInfo1.ref.spanningPairCount, 0);

  // Fragment supports both BP1 and BP2 of alt allele, but does not
  // support BP1 or BP2 of ref allele.
  SVSampleInfo sampleInfo2;
  fragmentEvidence.alt.bp1.isFragmentSupport = true;
  fragmentEvidence.alt.bp2.isFragmentSupport = true;
  addSpanningPairSupport(fragmentEvidence, sampleInfo2);
  BOOST_REQUIRE_EQUAL(sampleInfo2.alt.spanningPairCount, 1);
  BOOST_REQUIRE_EQUAL(sampleInfo2.ref.spanningPairCount, 0);

  // Fragment supports only BP1 of alt allele, but does not
  // support BP1 or BP2 of ref allele.
  SVSampleInfo sampleInfo3;
  fragmentEvidence.alt.bp1.isFragmentSupport = true;
  fragmentEvidence.alt.bp2.isFragmentSupport = false;
  addSpanningPairSupport(fragmentEvidence, sampleInfo3);
  BOOST_REQUIRE_EQUAL(sampleInfo3.alt.spanningPairCount, 1);
  BOOST_REQUIRE_EQUAL(sampleInfo3.ref.spanningPairCount, 0);

  // Fragment supports only BP2 of alt allele, but does not
  // support BP1 or BP2 of ref allele.
  SVSampleInfo sampleInfo4;
  fragmentEvidence.alt.bp1.isFragmentSupport = false;
  fragmentEvidence.alt.bp2.isFragmentSupport = true;
  addSpanningPairSupport(fragmentEvidence, sampleInfo4);
  BOOST_REQUIRE_EQUAL(sampleInfo4.alt.spanningPairCount, 1);
  BOOST_REQUIRE_EQUAL(sampleInfo4.ref.spanningPairCount, 0);

  // Fragment supports alt allele as well as supports both BP1 and BP2
  // of ref allele.
  SVSampleInfo sampleInfo5;
  fragmentEvidence.ref.bp1.isFragmentSupport = true;
  fragmentEvidence.ref.bp2.isFragmentSupport = true;
  addSpanningPairSupport(fragmentEvidence, sampleInfo5);
  BOOST_REQUIRE_EQUAL(sampleInfo5.alt.spanningPairCount, 1);
  BOOST_REQUIRE_EQUAL(sampleInfo5.ref.spanningPairCount, 1);

  // Fragment supports alt allele as well as supports only BP1
  // of ref allele.
  SVSampleInfo sampleInfo6;
  fragmentEvidence.ref.bp1.isFragmentSupport = true;
  fragmentEvidence.ref.bp2.isFragmentSupport = false;
  addSpanningPairSupport(fragmentEvidence, sampleInfo6);
  BOOST_REQUIRE_EQUAL(sampleInfo6.alt.spanningPairCount, 1);
  BOOST_REQUIRE_EQUAL(sampleInfo6.ref.spanningPairCount, 1);

  // Fragment supports alt allele as well as supports only BP2
  // of ref allele.
  SVSampleInfo sampleInfo7;
  fragmentEvidence.ref.bp1.isFragmentSupport = false;
  fragmentEvidence.ref.bp2.isFragmentSupport = true;
  addSpanningPairSupport(fragmentEvidence, sampleInfo7);
  BOOST_REQUIRE_EQUAL(sampleInfo7.alt.spanningPairCount, 1);
  BOOST_REQUIRE_EQUAL(sampleInfo7.ref.spanningPairCount, 1);
}

// Following cases need to be tested for confidentSemiMappedSpanningPairCount,
// confidentSpanningPairCount of alt allele and for confidentSemiMappedSpanningPairCount,
// confidentSpanningPairCount, confidentSplitReadAndPairCountRefBp1,
// confidentSplitReadAndPairCountRefBp2 of ref allele:
// 1. When a fragment does not support BP1 or BP2 of both ref and allele. So all the counts are 0.
// 2. If fragment size probability of alt allele is greater than fragment size probability of
//    ref allele and fraction of fragment size probability of alt allele greater than 0.9, then
//    confidentSemiMappedSpanningPairCount of alt allele is incremented by 1.
// 3. If case-2 condition is satisfied for fully mapped fragment where fully mapped fragment means
//    both read1 and read2 have a confident mapping, then
//    confidentSpanningPairCount of alt allele is incremented by 1.
// 4. If fragment size probability of ref allele is greater than fragment size probability of
//    alt allele and fraction of fragment size probability of ref allele greater than 0.9, then
//    ConfidentSemiMappedSpanningPairCount of ref allele is incremented by 1.
// 5. If case-4 condition is satisfied for fully mapped fragment where fully mapped fragment means
//    both read1 and read2 have a confident mapping, then
//    confidentSpanningPairCount of ref allele is incremented by 1.
// 6. If case-4 condition is satisfied for fully mapped fragment, then increment
// confidentSplitReadAndPairCountRefBp1
//    or confidentSplitReadAndPairCountRefBp2 or both by 1 based on whether this fragment supports this allele
//    on BP1 or BP2 or both respectively.
// 7. If fragment size probability = 0 for both alt and ref allele, it throws an exception.
BOOST_AUTO_TEST_CASE(test_addConservativeSpanningPairSupport)
{
  float fragLengthProb1(0.4);
  // Designed the case-1 where a fragment does not support any allele on BP1 or BP2
  SVFragmentEvidence fragmentEvidence1;
  SVSampleInfo       sampleInfo1;
  addConservativeSpanningPairSupport(fragmentEvidence1, sampleInfo1);
  BOOST_REQUIRE_EQUAL(sampleInfo1.alt.confidentSemiMappedSpanningPairCount, 0);
  BOOST_REQUIRE_EQUAL(sampleInfo1.alt.confidentSpanningPairCount, 0);
  BOOST_REQUIRE_EQUAL(sampleInfo1.ref.confidentSemiMappedSpanningPairCount, 0);
  BOOST_REQUIRE_EQUAL(sampleInfo1.ref.confidentSpanningPairCount, 0);
  BOOST_REQUIRE_EQUAL(sampleInfo1.ref.confidentSplitReadAndPairCountRefBp1, 0);
  BOOST_REQUIRE_EQUAL(sampleInfo1.ref.confidentSplitReadAndPairCountRefBp2, 0);

  // Designed the case-2
  // Here altFragProbability = 0.4 and refFragProbability = 0
  // so altFragProbability > refFragProbability
  // and altFragProbability/(altFragProbability + refFragProbability) = 1 > 0.9
  SVSampleInfo sampleInfo2;
  fragmentEvidence1.alt.bp1.isFragmentSupport = true;
  fragmentEvidence1.alt.bp1.fragLengthProb    = fragLengthProb1;
  addConservativeSpanningPairSupport(fragmentEvidence1, sampleInfo2);
  BOOST_REQUIRE_EQUAL(sampleInfo2.alt.confidentSemiMappedSpanningPairCount, 1);
  BOOST_REQUIRE_EQUAL(sampleInfo2.alt.confidentSpanningPairCount, 0);
  BOOST_REQUIRE_EQUAL(sampleInfo2.ref.confidentSemiMappedSpanningPairCount, 0);
  BOOST_REQUIRE_EQUAL(sampleInfo2.ref.confidentSpanningPairCount, 0);
  BOOST_REQUIRE_EQUAL(sampleInfo2.ref.confidentSplitReadAndPairCountRefBp1, 0);
  BOOST_REQUIRE_EQUAL(sampleInfo2.ref.confidentSplitReadAndPairCountRefBp2, 0);

  // Designed the case-3
  // Here altFragProbability = 0.4 and refFragProbability = 0
  // so altFragProbability > refFragProbability
  // and altFragProbability/(altFragProbability + refFragProbability) = 1 > 0.9
  // Also read1 and read2 both are anchored read means both have confident mapping.
  SVSampleInfo sampleInfo3;
  fragmentEvidence1.alt.bp1.isFragmentSupport = true;
  fragmentEvidence1.alt.bp1.fragLengthProb    = fragLengthProb1;
  fragmentEvidence1.read1.isScanned           = true;
  fragmentEvidence1.read1.setAnchored(true);
  fragmentEvidence1.read2.isScanned = true;
  fragmentEvidence1.read2.setAnchored(true);
  addConservativeSpanningPairSupport(fragmentEvidence1, sampleInfo3);
  BOOST_REQUIRE_EQUAL(sampleInfo3.alt.confidentSemiMappedSpanningPairCount, 1);
  BOOST_REQUIRE_EQUAL(sampleInfo3.alt.confidentSpanningPairCount, 1);
  BOOST_REQUIRE_EQUAL(sampleInfo3.ref.confidentSemiMappedSpanningPairCount, 0);
  BOOST_REQUIRE_EQUAL(sampleInfo3.ref.confidentSpanningPairCount, 0);
  BOOST_REQUIRE_EQUAL(sampleInfo3.ref.confidentSplitReadAndPairCountRefBp1, 0);
  BOOST_REQUIRE_EQUAL(sampleInfo3.ref.confidentSplitReadAndPairCountRefBp2, 0);

  // Designed the case-4
  // Here refFragProbability = 0.4 and altFragProbability = 0
  // so refFragProbability > altFragProbability
  // and refFragProbability/(altFragProbability + refFragProbability) = 1 > 0.9
  SVFragmentEvidence fragmentEvidence2;
  SVSampleInfo       sampleInfo4;
  fragmentEvidence2.ref.bp1.isFragmentSupport = true;
  fragmentEvidence2.ref.bp1.fragLengthProb    = fragLengthProb1;
  addConservativeSpanningPairSupport(fragmentEvidence2, sampleInfo4);
  BOOST_REQUIRE_EQUAL(sampleInfo4.alt.confidentSemiMappedSpanningPairCount, 0);
  BOOST_REQUIRE_EQUAL(sampleInfo4.alt.confidentSpanningPairCount, 0);
  BOOST_REQUIRE_EQUAL(sampleInfo4.ref.confidentSemiMappedSpanningPairCount, 1);
  BOOST_REQUIRE_EQUAL(sampleInfo4.ref.confidentSpanningPairCount, 0);
  BOOST_REQUIRE_EQUAL(sampleInfo4.ref.confidentSplitReadAndPairCountRefBp1, 0);
  BOOST_REQUIRE_EQUAL(sampleInfo4.ref.confidentSplitReadAndPairCountRefBp2, 0);

  // Designed the case-5
  // Here refFragProbability = 0.4 and altFragProbability = 0
  // so refFragProbability > altFragProbability
  // and refFragProbability/(altFragProbability + refFragProbability) = 1 > 0.9
  // Also read1 and read2 both are anchored read means both have confident mapping.
  SVSampleInfo sampleInfo5;
  fragmentEvidence2.ref.bp1.isFragmentSupport = true;
  fragmentEvidence2.ref.bp1.fragLengthProb    = fragLengthProb1;
  fragmentEvidence2.read1.isScanned           = true;
  fragmentEvidence2.read1.setAnchored(true);
  fragmentEvidence2.read2.isScanned = true;
  fragmentEvidence2.read2.setAnchored(true);
  addConservativeSpanningPairSupport(fragmentEvidence2, sampleInfo5);
  BOOST_REQUIRE_EQUAL(sampleInfo5.alt.confidentSemiMappedSpanningPairCount, 0);
  BOOST_REQUIRE_EQUAL(sampleInfo5.alt.confidentSpanningPairCount, 0);
  BOOST_REQUIRE_EQUAL(sampleInfo5.ref.confidentSemiMappedSpanningPairCount, 1);
  BOOST_REQUIRE_EQUAL(sampleInfo5.ref.confidentSpanningPairCount, 1);
  // Designed case-6 where fragment supports BP1 of ref allele
  BOOST_REQUIRE_EQUAL(sampleInfo5.ref.confidentSplitReadAndPairCountRefBp1, 1);
  BOOST_REQUIRE_EQUAL(sampleInfo5.ref.confidentSplitReadAndPairCountRefBp2, 0);

  // Designed the case-5
  // Here refFragProbability = 0.4 and altFragProbability = 0
  // so refFragProbability > altFragProbability
  // and refFragProbability/(altFragProbability + refFragProbability) = 1 > 0.9
  // Also read1 and read2 both are anchored read.
  // Designed case-6 where fragment supports BP2 of ref allele
  SVSampleInfo sampleInfo6;
  fragmentEvidence2.ref.bp1.isFragmentSupport = false;
  fragmentEvidence2.ref.bp2.isFragmentSupport = true;
  fragmentEvidence2.ref.bp2.fragLengthProb    = fragLengthProb1;
  fragmentEvidence2.read1.isScanned           = true;
  fragmentEvidence2.read1.setAnchored(true);
  fragmentEvidence2.read2.isScanned = true;
  fragmentEvidence2.read2.setAnchored(true);
  addConservativeSpanningPairSupport(fragmentEvidence2, sampleInfo6);
  BOOST_REQUIRE_EQUAL(sampleInfo6.alt.confidentSemiMappedSpanningPairCount, 0);
  BOOST_REQUIRE_EQUAL(sampleInfo6.alt.confidentSpanningPairCount, 0);
  BOOST_REQUIRE_EQUAL(sampleInfo6.ref.confidentSemiMappedSpanningPairCount, 1);
  BOOST_REQUIRE_EQUAL(sampleInfo6.ref.confidentSpanningPairCount, 1);
  BOOST_REQUIRE_EQUAL(sampleInfo6.ref.confidentSplitReadAndPairCountRefBp1, 0);
  // fragment supports BP2 of ref allele
  BOOST_REQUIRE_EQUAL(sampleInfo6.ref.confidentSplitReadAndPairCountRefBp2, 1);

  // Designed the case-7 where both alt and ref fragment size probability = 0.
  // It will throw an exception
  fragmentEvidence2.ref.bp2.fragLengthProb = 0;
  BOOST_CHECK_THROW(
      addConservativeSpanningPairSupport(fragmentEvidence2, sampleInfo6), illumina::common::GeneralException);
}

// For each fragment evidence api calculates the following stats and combines it:
// 1. confidentSpanningPairCount, confidentSpanningPairCount, confidentSemiMappedSpanningPairCount which
//    are explained in test_addConservativeSpanningPairSupport
// 2. confidentSplitReadCount which is explained in test_addConservativeSplitReadSupport
// 3. spanningPairCount which is explained in test_addSpanningPairSupport
// Test whether the api does aggregation of fragments evidences.
BOOST_AUTO_TEST_CASE(test_getSampleCounts)
{
  // sample evidence information
  SVEvidence::evidenceTrack_t samples;
  SVFragmentEvidence          fragmentEvidence1;
  SVSampleInfo                sampleInfo1;
  fragmentEvidence1.alt.bp1.isFragmentSupport = true;
  fragmentEvidence1.alt.bp1.fragLengthProb    = 0.4f;
  fragmentEvidence1.read1.isScanned           = true;
  fragmentEvidence1.read1.setAnchored(true);
  fragmentEvidence1.read2.isScanned = true;
  fragmentEvidence1.read2.setAnchored(true);
  samples["Fragment-1"] = fragmentEvidence1;
  // Here for fragment-1, altFragProbability = 0.4 and refFragProbability = 0
  // so altFragProbability > refFragProbability
  // and altFragProbability/(altFragProbability + refFragProbability) = 1 > 0.9.
  // So confidentSemiMappedSpanningPairCount is incremented by 1.
  // Also read1 and read2 both are anchored read. So,  confidentSpanningPairCount is
  // incremented by 1.
  // For details, you can see test_addConservativeSpanningPairSupport.
  getSampleCounts(samples, sampleInfo1);
  BOOST_REQUIRE_EQUAL(sampleInfo1.alt.confidentSpanningPairCount, 1);
  BOOST_REQUIRE_EQUAL(sampleInfo1.alt.confidentSemiMappedSpanningPairCount, 1);
  SVSampleInfo       sampleInfo2;
  SVFragmentEvidence fragmentEvidence2;
  fragmentEvidence2.alt.bp1.read1.isSplitSupport = true;
  fragmentEvidence2.alt.bp1.read1.splitLnLhood   = -7.9;
  fragmentEvidence2.alt.bp2.read1.splitLnLhood   = -8.9;
  fragmentEvidence2.ref.bp1.read1.splitLnLhood   = -17.2;
  fragmentEvidence2.ref.bp2.read1.splitLnLhood   = -18.9;
  samples["Fragment-2"]                          = fragmentEvidence2;
  // For fragment-2, altLnLhood = max(bp1 splitLnLhood of alt, bp2 splitLnLhood of alt) = -7.9
  // and refLnLhood = max(bp1 splitLnLhood of ref, bp2 splitLnLhood of ref) = -17.2
  // Here, altLnLhood > refLnLhood. So normalized probability(based on test_lnToProb)
  // of altLnLhood  = 1/exp(refLnLhood - altLnLhood)
  //                = 0.999909 which is greater than 0.999f
  // So it is a confident split read evidence for alt allele.
  // Detailed has been explained in test_addConservativeSplitReadSupport
  getSampleCounts(samples, sampleInfo2);
  // following two cases are coming from fragment-1 in test_addConservativeSpanningPairSupport
  BOOST_REQUIRE_EQUAL(sampleInfo2.alt.confidentSpanningPairCount, 1);
  BOOST_REQUIRE_EQUAL(sampleInfo2.alt.confidentSemiMappedSpanningPairCount, 1);
  // This one is coming from fragment-2 in test_addConservativeSplitReadSupport
  BOOST_REQUIRE_EQUAL(sampleInfo2.alt.confidentSplitReadCount, 1);

  SVSampleInfo       sampleInfo3;
  SVFragmentEvidence fragmentEvidence3;
  fragmentEvidence3.alt.bp1.isFragmentSupport = true;
  fragmentEvidence3.alt.bp1.fragLengthProb    = 0.6f;
  fragmentEvidence3.ref.bp1.fragLengthProb    = 0.04f;
  samples["Fragment-3"]                       = fragmentEvidence3;
  // Here for fragment-3, altFragProbability = 0.6 and refFragProbability = 0.04
  // so altFragProbability > refFragProbability
  // and altFragProbability/(altFragProbability + refFragProbability) = 0.935 > 0.9.
  // So confidentSemiMappedSpanningPairCount is incremented by 1.
  getSampleCounts(samples, sampleInfo3);
  // So total number of confidentSemiMappedSpanningPairCount is 2 (fragment-1 and fragment-3)
  BOOST_REQUIRE_EQUAL(sampleInfo3.alt.confidentSemiMappedSpanningPairCount, 2);
  // total number of confidentSpanningPairCount is 1 (fragment-1)
  BOOST_REQUIRE_EQUAL(sampleInfo3.alt.confidentSpanningPairCount, 1);
  // total number of confidentSplitReadCount is 1 (fragment-2)
  BOOST_REQUIRE_EQUAL(sampleInfo3.alt.confidentSplitReadCount, 1);
  // total number of spanning pair count is 1 (fragment-1 and fragment-3) as it supports alt allele on BP1.
  BOOST_REQUIRE_EQUAL(sampleInfo3.alt.spanningPairCount, 2);
}

// For each sample api calculates the following stats:
// 1. confidentSpanningPairCount, confidentSpanningPairCount, confidentSemiMappedSpanningPairCount which
//    are explained in test_addConservativeSpanningPairSupport
// 2. confidentSplitReadCount which is explained in test_addConservativeSplitReadSupport
// 3. spanningPairCount which is explained in test_addSpanningPairSupport
// Test whether the api does stats computation  of fragments evidences of each sample.
BOOST_AUTO_TEST_CASE(test_getSVSupportSummary)
{
  // Data for sample-1
  SVEvidence::evidenceTrack_t sample1;
  // Here for fragment-1, altFragProbability = 0.4 and refFragProbability = 0
  // so altFragProbability > refFragProbability
  // and altFragProbability/(altFragProbability + refFragProbability) = 1 > 0.9.
  // So confidentSemiMappedSpanningPairCount is incremented by 1.
  // Also read1 and read2 both are anchored read. So,  confidentSpanningPairCount is
  // incremented by 1.
  // For details, you can see test_addConservativeSpanningPairSupport.
  SVFragmentEvidence fragmentEvidence1;
  fragmentEvidence1.alt.bp1.isFragmentSupport = true;
  fragmentEvidence1.alt.bp1.fragLengthProb    = 0.4f;
  fragmentEvidence1.read1.isScanned           = true;
  fragmentEvidence1.read1.setAnchored(true);
  fragmentEvidence1.read2.isScanned = true;
  fragmentEvidence1.read2.setAnchored(true);
  sample1["Fragment-1"] = fragmentEvidence1;
  // For fragment-2, altLnLhood = max(bp1 splitLnLhood of alt, bp2 splitLnLhood of alt) = -7.9
  // and refLnLhood = max(bp1 splitLnLhood of ref, bp2 splitLnLhood of ref) = -17.2
  // Here, altLnLhood > refLnLhood. So normalized probability(based on test_lnToProb)
  // of altLnLhood  = 1/exp(refLnLhood - altLnLhood)
  //                = 0.999909 which is greater than 0.999f
  // So it is a confident split read evidence for alt allele.
  // Detailed has been explained in test_addConservativeSplitReadSupport
  SVFragmentEvidence fragmentEvidence2;
  fragmentEvidence2.alt.bp1.read1.isSplitSupport = true;
  fragmentEvidence2.alt.bp1.read1.splitLnLhood   = -7.9;
  fragmentEvidence2.alt.bp2.read1.splitLnLhood   = -8.9;
  fragmentEvidence2.ref.bp1.read1.splitLnLhood   = -17.2;
  fragmentEvidence2.ref.bp2.read1.splitLnLhood   = -18.9;
  sample1["Fragment-2"]                          = fragmentEvidence2;
  SVEvidence evidence;
  evidence.samples.push_back(sample1);

  // Data for sample-2
  SVEvidence::evidenceTrack_t sample2;
  // Here for fragment-1 of sample-2, altFragProbability = 0.6 and refFragProbability = 0.04
  // so altFragProbability > refFragProbability
  // and altFragProbability/(altFragProbability + refFragProbability) = 0.935 > 0.9.
  // So confidentSemiMappedSpanningPairCount is incremented by 1.
  SVFragmentEvidence fragmentEvidence3;
  fragmentEvidence3.alt.bp1.isFragmentSupport = true;
  fragmentEvidence3.alt.bp1.fragLengthProb    = 0.6f;
  fragmentEvidence3.ref.bp1.fragLengthProb    = 0.04f;
  sample2["Fragment-1"]                       = fragmentEvidence3;
  // For fragment-2 of sample-2, refLnLhood = max(bp1 splitLnLhood of ref, bp2 splitLnLhood of ref) = -7.9
  // and altLnLhood = max(bp1 splitLnLhood of alt, bp2 splitLnLhood of alt) = -17.2
  // Here, refLnLhood > altLnLhood. So normalized probability(based on test_lnToProb)
  // of refLnLhood  = 1/exp(altLnLhood - refLnLhood)
  //                = 0.999909 which is greater than 0.999f
  // So it is a confident split read evidence for ref allele.
  // Detailed has been explained in test_addConservativeSplitReadSupport
  SVFragmentEvidence fragmentEvidence4;
  fragmentEvidence4.ref.bp1.read1.isSplitSupport = true;
  fragmentEvidence4.ref.bp1.read1.splitLnLhood   = -7.9;
  fragmentEvidence4.ref.bp2.read1.splitLnLhood   = -8.9;
  fragmentEvidence4.alt.bp1.read1.splitLnLhood   = -17.2;
  fragmentEvidence4.alt.bp2.read1.splitLnLhood   = -18.9;
  sample2["Fragment-2"]                          = fragmentEvidence4;
  evidence.samples.push_back(sample2);
  SVScoreInfo scoreInfo;
  scoreInfo.samples.resize(2);
  getSVSupportSummary(evidence, scoreInfo);

  // check informations for sample-1
  // Total number of confidentSemiMappedSpanningPairCount and confidentSpanningPairCount
  // of alt allele is 1 (for fragment-1).
  BOOST_REQUIRE_EQUAL(scoreInfo.samples[0].alt.confidentSpanningPairCount, 1);
  BOOST_REQUIRE_EQUAL(scoreInfo.samples[0].alt.confidentSemiMappedSpanningPairCount, 1);
  // Total number of confidentSplitReadCount of alt allele is 1 (for fragment-2).
  BOOST_REQUIRE_EQUAL(scoreInfo.samples[0].alt.confidentSplitReadCount, 1);
  // check informations for sample-2
  // Total number of confidentSemiMappedSpanningPairCount alt allele is 1 (for fragment-1).
  BOOST_REQUIRE_EQUAL(scoreInfo.samples[1].alt.confidentSemiMappedSpanningPairCount, 1);
  BOOST_REQUIRE_EQUAL(scoreInfo.samples[1].alt.confidentSpanningPairCount, 0);
  // Total number of confidentSplitReadCount of ref allele is 1 (for fragment-2).
  BOOST_REQUIRE_EQUAL(scoreInfo.samples[1].ref.confidentSplitReadCount, 1);
}

// If there is a conflict of split support and fragment support for a read-pair then
// priority should be given to split evidence over fragment support in some scenarioes. In that case, api
// clears the fragment support evidence. Let's consider AFP, RFP, ASL, RSL are alt allele fragment size
// probability, ref allele fragment size probability, altSplitLnLhood and refSplitLnLhood respectively, then
// test
// whether api clears fragment support evidence for the following cases:
// 1. When RSL > ASL at read1, but AFP > RFP
// 2. When ASL > RSL at read1, but RFP > AFP
// 3. When RSL > ASL at read2, but AFP > RFP
// 4. When ASL > RSL at read2, but RFP > AFP
BOOST_AUTO_TEST_CASE(test_resolvePairSplitConflicts)
{
  float       fragLengthProb1(0.4);
  SVCandidate candidate;
  candidate.setPrecise();
  SVEvidence::evidenceTrack_t evidenceTrack;
  // Designed the case-1 where read-1 supports split evidence
  std::string        fragLabel = "frag-1";
  SVFragmentEvidence fragmentEvidence1;
  // Here for this fragment, altFragProbability = 0.4 and refFragProbability = 0
  // so altFragProbability > refFragProbability
  fragmentEvidence1.alt.bp1.isFragmentSupport    = true;
  fragmentEvidence1.alt.bp1.fragLengthProb       = fragLengthProb1;
  fragmentEvidence1.ref.bp1.read1.isSplitSupport = true;
  // For this fragment , refLnLhood = max(bp1 splitLnLhood of ref, bp2 splitLnLhood of ref) = -7.9
  // and altLnLhood = max(bp1 splitLnLhood of alt, bp2 splitLnLhood of alt) = -17.2
  // Here, refLnLhood > altLnLhood.
  fragmentEvidence1.ref.bp1.read1.splitLnLhood = -7.9;
  fragmentEvidence1.ref.bp2.read1.splitLnLhood = -8.9;
  fragmentEvidence1.alt.bp1.read1.splitLnLhood = -17.2;
  fragmentEvidence1.alt.bp2.read1.splitLnLhood = -18.9;
  evidenceTrack[fragLabel]                     = fragmentEvidence1;
  SVEvidence evidence;
  evidence.samples.resize(1);
  evidence.samples.push_back(evidenceTrack);
  // So RSL > ASL at read1, but AFP > RFP. Clear the fragment support evidence.
  resolvePairSplitConflicts(candidate, evidence);
  BOOST_REQUIRE(!evidence.getSampleEvidence(0)[fragLabel].alt.bp1.isFragmentSupport);
  BOOST_REQUIRE_EQUAL(evidence.getSampleEvidence(0)[fragLabel].alt.bp1.fragLengthProb, 0);

  // Designed the case-2 where read-1 supports split evidence
  SVFragmentEvidence fragmentEvidence2;
  // Here for this fragment, refFragProbability = 0.4 and altFragProbability = 0
  // so refFragProbability > altFragProbability
  fragmentEvidence2.ref.bp1.isFragmentSupport    = true;
  fragmentEvidence2.ref.bp1.fragLengthProb       = fragLengthProb1;
  fragmentEvidence2.alt.bp1.read1.isSplitSupport = true;
  // For this fragment , altLnLhood = max(bp1 splitLnLhood of alt, bp2 splitLnLhood of alt) = -7.9
  // and refLnLhood = max(bp1 splitLnLhood of ref, bp2 splitLnLhood of ref) = -17.2
  // Here, altLnLhood > refLnLhood.
  fragmentEvidence2.alt.bp1.read1.splitLnLhood = -7.9;
  fragmentEvidence2.alt.bp2.read1.splitLnLhood = -8.9;
  fragmentEvidence2.ref.bp1.read1.splitLnLhood = -17.2;
  fragmentEvidence2.ref.bp2.read1.splitLnLhood = -18.9;
  evidenceTrack[fragLabel]                     = fragmentEvidence2;
  evidence.samples.clear();
  evidence.samples.push_back(evidenceTrack);
  // So ASL > RSL at read1, but RFP > AFP. Clear the fragment support evidence.
  resolvePairSplitConflicts(candidate, evidence);
  BOOST_REQUIRE(!evidence.getSampleEvidence(0)[fragLabel].ref.bp1.isFragmentSupport);
  BOOST_REQUIRE_EQUAL(evidence.getSampleEvidence(0)[fragLabel].ref.bp1.fragLengthProb, 0);

  // Designed the case-3 where read-2 supports split evidence
  SVFragmentEvidence fragmentEvidence3;
  // Here for this fragment, altFragProbability = 0.4 and refFragProbability = 0
  // so altFragProbability > refFragProbability
  fragmentEvidence3.alt.bp1.isFragmentSupport = true;
  fragmentEvidence3.alt.bp1.fragLengthProb    = fragLengthProb1;
  // For this fragment , refLnLhood = max(bp1 splitLnLhood of ref, bp2 splitLnLhood of ref) = -7.9
  // and altLnLhood = max(bp1 splitLnLhood of alt, bp2 splitLnLhood of alt) = -17.2
  // Here, refLnLhood > altLnLhood.
  fragmentEvidence3.alt.bp1.read2.isSplitSupport = true;
  fragmentEvidence3.ref.bp1.read1.splitLnLhood   = -7.9;
  fragmentEvidence3.ref.bp2.read1.splitLnLhood   = -8.9;
  fragmentEvidence3.alt.bp1.read1.splitLnLhood   = -17.2;
  fragmentEvidence1.alt.bp2.read1.splitLnLhood   = -18.9;
  evidenceTrack[fragLabel]                       = fragmentEvidence3;
  evidence.samples.clear();
  evidence.samples.push_back(evidenceTrack);
  // So RSL > ASL at read1, but AFP > RFP. Clear the fragment support evidence.
  resolvePairSplitConflicts(candidate, evidence);
  BOOST_REQUIRE(!evidence.getSampleEvidence(0)[fragLabel].ref.bp1.isFragmentSupport);
  BOOST_REQUIRE_EQUAL(evidence.getSampleEvidence(0)[fragLabel].ref.bp1.fragLengthProb, 0);

  // Designed the case-4 where read-2 supports split evidence
  SVFragmentEvidence fragmentEvidence4;
  // Here for this fragment, refFragProbability = 0.4 and altFragProbability = 0
  // so refFragProbability > altFragProbability
  fragmentEvidence4.ref.bp1.isFragmentSupport    = true;
  fragmentEvidence4.ref.bp1.fragLengthProb       = fragLengthProb1;
  fragmentEvidence4.ref.bp1.read2.isSplitSupport = true;
  // For this tragment , altLnLhood = max(bp1 splitLnLhood of alt, bp2 splitLnLhood of alt) = -7.9
  // and refLnLhood = max(bp1 splitLnLhood of ref, bp2 splitLnLhood of ref) = -17.2
  // Here, altLnLhood > refLnLhood.
  fragmentEvidence4.alt.bp1.read1.splitLnLhood = -7.9;
  fragmentEvidence4.alt.bp2.read1.splitLnLhood = -8.9;
  fragmentEvidence4.ref.bp1.read1.splitLnLhood = -17.2;
  fragmentEvidence4.ref.bp2.read1.splitLnLhood = -18.9;
  evidenceTrack[fragLabel]                     = fragmentEvidence4;
  evidence.samples.clear();
  evidence.samples.push_back(evidenceTrack);
  // So ASL > RSL at read1, but RFP > AFP. Clear the fragment support evidence.
  resolvePairSplitConflicts(candidate, evidence);
  BOOST_REQUIRE(!evidence.getSampleEvidence(0)[fragLabel].ref.bp1.isFragmentSupport);
  BOOST_REQUIRE_EQUAL(evidence.getSampleEvidence(0)[fragLabel].ref.bp1.fragLengthProb, 0);
}

// Test the following utility formula:
// log(selfChimeraProb.comp*fragProb + otherChimeraProb.prob)*power,
// where fragProb is
// 1. fragment-size probability of breakpoint-1, if a fragment supports only breakpoint-1
// 2. fragment-size probability of breakpoint-2, if a fragment supports only breakpoint-2
// 3. maximum of fragment-size probability of breakpoint-1 and breakpoint-2, if a fragment
//    supports both breakpoint-1 and breakpoint-2
// Here the value of power is 2.
BOOST_AUTO_TEST_CASE(test_incrementSpanningPairAlleleLnLhood)
{
  float               fragLengthProb1(0.4);
  float               fragLengthProb2(0.6);
  ProbSet             selfChimeraProb(0.5);
  ProbSet             otherChimeraProb(0.2);
  static const double eps(0.00000001);
  // Case-1 when fragment supports only BP1
  SVFragmentEvidenceAllele fragmentEvidenceAllele;
  fragmentEvidenceAllele.bp1.isFragmentSupport = true;
  fragmentEvidenceAllele.bp1.fragLengthProb    = fragLengthProb1;
  double bpLnLhood1Expected(-1.832581448847149);
  double bpLnLhood1(0);
  incrementSpanningPairAlleleLnLhood(
      selfChimeraProb, otherChimeraProb, fragmentEvidenceAllele, 2, bpLnLhood1);
  BOOST_REQUIRE_CLOSE(bpLnLhood1, bpLnLhood1Expected, eps);

  // Case-2 when fragment supports only BP2
  double bpLnLhood2Expected(-1.3862943134361756);
  double bpLnLhood2(0);
  fragmentEvidenceAllele.bp1.isFragmentSupport = false;
  fragmentEvidenceAllele.bp2.isFragmentSupport = true;
  fragmentEvidenceAllele.bp2.fragLengthProb    = fragLengthProb2;
  incrementSpanningPairAlleleLnLhood(
      selfChimeraProb, otherChimeraProb, fragmentEvidenceAllele, 2, bpLnLhood2);
  BOOST_REQUIRE_CLOSE(bpLnLhood2, bpLnLhood2Expected, eps);

  // Case-3 when fragment supports both BP1 and BP2
  double bpLnLhood3(0);
  fragmentEvidenceAllele.bp1.isFragmentSupport = true;
  fragmentEvidenceAllele.bp2.isFragmentSupport = true;
  fragmentEvidenceAllele.bp2.fragLengthProb    = fragLengthProb2;
  incrementSpanningPairAlleleLnLhood(
      selfChimeraProb, otherChimeraProb, fragmentEvidenceAllele, 2, bpLnLhood3);
  BOOST_REQUIRE_CLOSE(bpLnLhood3, bpLnLhood2Expected, eps);
}

// fragLnLHood = log_sum(selfMapProb.lnComp+alignLnLhood), otherMapProb.lnProb),
//     log_sum(x1,x2) = x1 + log(1+exp(x2-x1)) if x2>x1
//                    = x2 + log(1+exp(x1-x2)) if x1>x2
//     lnComp = log(1 - lnProb)
// Following cases need to tested
// 1. When a read does not support split evidence on any of the breakpoints
// 2. When a read supports split evidence on BP1 but does not support split evidence on BP2
// 3. When a read supports split evidence on BP2 but does not support split evidence on BP1
// 4. When a read supports split evidence on both the breakpoints and splitLnLhood of BP1
//    greater than splitLnLhood of BP2
// 5. When a read supports split evidence on both the breakpoints and splitLnLhood of BP2
//    greater than splitLnLhood of BP1
BOOST_AUTO_TEST_CASE(test_incrementAlleleSplitReadLhood)
{
  ProbSet             selfChimeraProb(0.5);
  ProbSet             otherChimeraProb(0.2);
  static const double eps(0.00000001);
  // Designed the case-1 where read-1 does not support any of the breakpoints.
  // So alignLnLhood = splitLnLhood of BP2 = -5
  SVFragmentEvidenceAllele fragmentEvidenceAllele;
  fragmentEvidenceAllele.bp1.read1.splitLnLhood = -6;
  fragmentEvidenceAllele.bp2.read1.splitLnLhood = -5;
  bool   isReadEvaluated;
  double expectedValue1(-1.5927333463365996);
  double expectedValue2(-1.6032601537000746);
  BOOST_REQUIRE_CLOSE(
      incrementAlleleSplitReadLhood(
          selfChimeraProb,
          otherChimeraProb,
          fragmentEvidenceAllele,
          0.5,
          fragmentEvidenceAllele.isAnySplitReadSupport(true),
          true,
          isReadEvaluated),
      expectedValue1,
      eps);

  // Designed the case-2 where read-1 supports split evidence on BP1 and does not support split evidence
  // on BP2.
  // So alignLnLhood = splitLnLhood of BP1 = -6
  fragmentEvidenceAllele.bp1.read1.isSplitSupport = true;
  BOOST_REQUIRE_CLOSE(
      incrementAlleleSplitReadLhood(
          selfChimeraProb,
          otherChimeraProb,
          fragmentEvidenceAllele,
          0.5,
          fragmentEvidenceAllele.isAnySplitReadSupport(true),
          true,
          isReadEvaluated),
      expectedValue2,
      eps);

  // Designed the case-3 where read-1 supports split evidence on BP2 and does not support split evidence
  // on BP1.
  // So alignLnLhood = splitLnLhood of BP2 = -5
  fragmentEvidenceAllele.bp1.read1.isSplitSupport = false;
  fragmentEvidenceAllele.bp2.read1.isSplitSupport = true;
  BOOST_REQUIRE_CLOSE(
      incrementAlleleSplitReadLhood(
          selfChimeraProb,
          otherChimeraProb,
          fragmentEvidenceAllele,
          0.5,
          fragmentEvidenceAllele.isAnySplitReadSupport(true),
          true,
          isReadEvaluated),
      expectedValue1,
      eps);

  // Designed the case-4 when read-1 supports split evidence on both the breakpoints and splitLnLhood of BP1
  // greater than splitLnLhood of BP2. Here alignLnLhood = -5
  fragmentEvidenceAllele.bp1.read1.isSplitSupport = true;
  fragmentEvidenceAllele.bp2.read1.isSplitSupport = true;
  fragmentEvidenceAllele.bp1.read1.splitLnLhood   = -5;
  fragmentEvidenceAllele.bp2.read1.splitLnLhood   = -6;
  BOOST_REQUIRE_CLOSE(
      incrementAlleleSplitReadLhood(
          selfChimeraProb,
          otherChimeraProb,
          fragmentEvidenceAllele,
          0.5,
          fragmentEvidenceAllele.isAnySplitReadSupport(true),
          true,
          isReadEvaluated),
      expectedValue1,
      eps);
  // Designed case-5 when read-1 supports split evidence on both the breakpoints and splitLnLhood of BP2
  // greater than splitLnLhood of BP1. Here alignLnLhood = -5
  fragmentEvidenceAllele.bp1.read1.splitLnLhood = -6;
  fragmentEvidenceAllele.bp2.read1.splitLnLhood = -5;
  BOOST_REQUIRE_CLOSE(
      incrementAlleleSplitReadLhood(
          selfChimeraProb,
          otherChimeraProb,
          fragmentEvidenceAllele,
          0.5,
          fragmentEvidenceAllele.isAnySplitReadSupport(true),
          true,
          isReadEvaluated),
      expectedValue1,
      eps);
}

// Following cases need to be tested
// 1. When a read does not support split evidence on any of the breakpoints
// 2. When a read does not support any tier2 split evidence on any of the breakpoints
// 3. When a read supports split evidence on one of the breakpoints
BOOST_AUTO_TEST_CASE(test_incrementSplitReadLhood)
{
  ProbSet             refMapProb(0.5);
  ProbSet             altMapProb(0.2);
  static const double eps(0.00000001);

  // Designed the case-1. In this case read-1 does not support split evidence on any of the breakpoints.
  SVFragmentEvidence fragmentEvidence;
  fragmentEvidence.alt.bp1.read1.splitLnLhood = -6;
  fragmentEvidence.alt.bp2.read1.splitLnLhood = -5;
  bool        isReadEvaluated;
  std::string fragLabel = "Frag-1";
  double      refSplitLnLhood(-10.3);
  double      altSplitLnLhood(-20.4);
  incrementSplitReadLhood(
      fragLabel,
      fragmentEvidence,
      refMapProb,
      altMapProb,
      false,
      true,
      refSplitLnLhood,
      altSplitLnLhood,
      isReadEvaluated);
  BOOST_REQUIRE_CLOSE(refSplitLnLhood, -10.3, eps);
  BOOST_REQUIRE_CLOSE(altSplitLnLhood, -20.4, eps);

  // Designed the case-2 when In this case read-1 does not support any tier2 split evidence
  // on any of the breakpoints.
  incrementSplitReadLhood(
      fragLabel,
      fragmentEvidence,
      refMapProb,
      altMapProb,
      true,
      true,
      refSplitLnLhood,
      altSplitLnLhood,
      isReadEvaluated);
  BOOST_REQUIRE_CLOSE(refSplitLnLhood, -10.3, eps);
  BOOST_REQUIRE_CLOSE(altSplitLnLhood, -20.4, eps);

  // Designed the case-3 read-1 supports split evidence on both the breakpoints and splitLnLhood of BP2
  // greater than splitLnLhood of BP1. Here alignLnLhood = -5.
  // So, fragLnLHood = log_sum(selfMapProb.lnComp+alignLnLhood), otherMapProb.lnProb),
  //     log_sum(x1,x2) = x1 + log(1+exp(x2-x1)) if x2>x1
  //                    = x2 + log(1+exp(x1-x2)) if x1>x2
  //     lnComp = log(1 - lnProb)
  // Detail has been explained in test_incrementAlleleSplitReadLhood
  fragmentEvidence.alt.bp1.read1.isSplitSupport = true;
  fragmentEvidence.alt.bp2.read1.isSplitSupport = true;
  incrementSplitReadLhood(
      fragLabel,
      fragmentEvidence,
      refMapProb,
      altMapProb,
      false,
      true,
      refSplitLnLhood,
      altSplitLnLhood,
      isReadEvaluated);
  BOOST_REQUIRE_CLOSE(refSplitLnLhood, -10.656674943938732, eps);
  BOOST_REQUIRE_CLOSE(altSplitLnLhood, -21.082424162960997, eps);
}

// Test the following cases:
// 1. When evaluation is done only for read-1. Then api should return sum of sum of the split
//    likelihood of read-1 and the fragment-size likelihood.
// 2. When evaluation is done only for read-2. Then api should return sum of the split
//    likelihood of read-2 and the fragment-size likelihood.
// 3. When evaluation is done based on both read-1 and read-2. Then api should return sum of
//    (maximum value of split likelihood of read-1 and
//    split likelihood of read-2) and the fragment-size likelihood.
// 4. When evaluation is done without read-1 and read-2. Then api should return value of the
//    fragment-size likelihood.
BOOST_AUTO_TEST_CASE(test_getFragLnLhood)
{
  // book keeping of lnLhood of read-1 and read-2
  AlleleLnLhood lnLhood;
  lnLhood.read1Split = -0.5;
  lnLhood.read2Split = -0.2;
  lnLhood.fragPair   = -0.4;
  static const double eps(0.00000001);

  // Designed the case-1 where evaluation is done based on read-1.
  BOOST_REQUIRE_CLOSE(getFragLnLhood(lnLhood, true, false), lnLhood.read1Split + lnLhood.fragPair, eps);

  // Designed the case-2 where evaluation is done based on read-2.
  BOOST_REQUIRE_CLOSE(getFragLnLhood(lnLhood, false, true), lnLhood.read2Split + lnLhood.fragPair, eps);

  // Designed the case-3 where evaluation is done based on both read-1 and read-2.
  BOOST_REQUIRE_CLOSE(getFragLnLhood(lnLhood, true, true), lnLhood.read2Split + lnLhood.fragPair, eps);

  // Designed the case-4 where evaluation is done without read-1 and read-2.
  BOOST_REQUIRE_CLOSE(getFragLnLhood(lnLhood, false, false), lnLhood.fragPair, eps);
}

// Spanning weight is calculated as,
// spanning weight = min(1, max(0, (val-_min)*_factor)),
// For large insertions(length > 100),
//      val = insert sequence size, _min = 100, _max = 150, & insertSizeRamp applied
// For deletion and small insertion,
//      val = center size of sv, _min = 300, _max = 500, & svSizeRamp applied
// factor = 1/(_max-_min),
// Following cases need to tested:
// 1. When SV type is not indel. In this case it will return 1.
// 2. When SV has only deletion
// 3. When SV has insertion of length less than 100
// 4  when SV has large insertion (length greater than 100)
BOOST_AUTO_TEST_CASE(test_getSpanningPairWeight)
{
  static const double eps(0.00000001);

  SVCandidate candidate;
  // Designed the case-1 Where candidate is not indel type SV.
  // Spanning weight is 1.
  candidate.bp1.interval = GenomeInterval(0, 1000, 1200);
  candidate.bp2.interval = GenomeInterval(1, 1000, 1200);
  BOOST_REQUIRE_CLOSE(getSpanningPairWeight(candidate), 1, eps);

  // Designed the case-2 where SV type is a deletion
  // _min = 300, _max = 500
  // val = center size = 550
  // _factor = 1/200 = 0.005
  // spanning weight = min(1, max(0, (val-_min)*_factor))
  //                 = min(1, max(0, (550-300)*0.005) = 1
  candidate.bp1.interval = GenomeInterval(0, 1000, 1200);
  candidate.bp2.interval = GenomeInterval(0, 500, 600);
  candidate.bp1.state    = SVBreakendState::LEFT_OPEN;
  candidate.bp2.state    = SVBreakendState::RIGHT_OPEN;
  BOOST_REQUIRE_CLOSE(getSpanningPairWeight(candidate), 1, eps);

  // Designed the case-3 Where SV has small insertion (size = 98)
  // _min = 300, _max = 500
  // val = center size = 115
  // _factor = 1/200 = 0.005
  // spanning weight = min(1, max(0, (val-_min)*_factor))
  //                 = min(1, max(0, (115-300)*0.005) = 0
  candidate.bp1.interval = GenomeInterval(0, 1000, 1200);
  candidate.bp2.interval = GenomeInterval(0, 980, 990);
  candidate.insertSeq =
      "GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGT"
      "ATTTTCGTCTGGGGGGTGTGCACGCGATAGCATTGCGAGACGC";
  BOOST_REQUIRE_CLOSE(getSpanningPairWeight(candidate), 0, eps);

  // Designed the case-4 Where SV has large insertion (size = 102)
  // _min = 100, _max = 150
  // val = insert sequence size = 102
  // _factor = 1/50 = 0.02
  // spanning weight = min(1, max(0, (val-_min)*_factor))
  //                 = min(1, max(0, (102-100)*0.02) = 0.04
  candidate.bp1.interval = GenomeInterval(0, 1000, 1200);
  candidate.bp2.interval = GenomeInterval(0, 980, 990);
  candidate.insertSeq =
      "GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGT"
      "ATTTTCGTCTGGGGGGTGTGCACGCGATAGCATTGCGAGACGCTGGA";
  BOOST_REQUIRE_CLOSE(getSpanningPairWeight(candidate), 0.04f, eps);
}

// Spanning weight is calculated as,
// spanning weight = min(1, max(0, (val-_min)*_factor)),
// For large SV,
//      val = center size of SV, _min = 5000, _max = 10000, & svSizeRamp applied
// factor = 1/(_max-_min),
// For interchrosomal SV event it will return 1.
BOOST_AUTO_TEST_CASE(test_largeNoiseSVPriorWeight)
{
  static const double eps(0.00000001);

  SVCandidate candidate;
  // Interchrosomal SV. So largeNoiseSVPriorWeight = 1.
  candidate.bp1.interval = GenomeInterval(0, 1000, 1200);
  candidate.bp2.interval = GenomeInterval(1, 1000, 1200);
  BOOST_REQUIRE_CLOSE(largeNoiseSVPriorWeight(candidate), 1, eps);

  // _min = 5000, _max = 10000
  // val = center size of candidate = 5925
  // _factor = 1/5000 = 0.0002
  // spanning weight = min(1, max(0, (val-_min)*_factor))
  //                 = min(1, max(0, (5925-5000)*0.02) = 0.185000002
  candidate.bp1.interval = GenomeInterval(0, 1000, 1200);
  candidate.bp2.interval = GenomeInterval(0, 7000, 7050);
  BOOST_REQUIRE_CLOSE(largeNoiseSVPriorWeight(candidate), 0.185000002f, eps);
}

// Test the following cases
// 1. API returns true if either of the following three conditions satisfied:
//       * If fragment supports ref or alt allele on any of the breakpoints
//       * If read-1 supports split evidence at any of the breakpoints
//       * If read-2 supports split evidence at any of the breakpoints
// 2. RefLnLhood and altLnLhood of a fragment evidence are calculated as follows:
//           RefLnLhood = log(selfChimeraProb.comp*fragSizeProb + otherChimeraProb.prob)*power,
//           where selfChimeraProb = 0.2 and otherChimeraProb = 0.3
//           AltLnLhood = log(selfChimeraProb.comp*fragSizeProb + otherChimeraProb.prob)*power,
//           where selfChimeraProb = 0.3 and otherChimeraProb = 0.2
//           power = spanningPairWeight, if both read-1 and read-2 have confident mapping wrt fragment support
//                 = spanningPairWeight * semiMappedPower, if fragment is semmi-mapped and fragment-size
//                   probability of alt allele is greater than ref allele.
//                 = 0, if fragment is semmi-mapped and fragment-size
//                   probability of ref allele is greater than alt allele.
BOOST_AUTO_TEST_CASE(test_getRefAltFromFrag)
{
  static const double eps(0.00000001);
  float               spanningPairWeight(0.2);
  double              semiMappedPower(2);
  ProbSet             refChimeraProb(0.2);
  ProbSet             altChimeraProb(0.3);
  ProbSet             refSplitMapProb(0.1);
  ProbSet             altSplitMapProb(0.5);
  std::string         fragLabel = "frag-1";
  SVFragmentEvidence  evidence;
  AlleleLnLhood       refLnLhoodSet;
  AlleleLnLhood       altLnLhoodSet;
  bool                isRead1Evaluated;
  bool                isRead2Evaluated;
  // Fragment does not support any of the allele on any of the breakpoints
  BOOST_REQUIRE(!getRefAltFromFrag(
      spanningPairWeight,
      semiMappedPower,
      refChimeraProb,
      altChimeraProb,
      refSplitMapProb,
      altSplitMapProb,
      true,
      fragLabel,
      evidence,
      refLnLhoodSet,
      altLnLhoodSet,
      isRead1Evaluated,
      isRead2Evaluated));

  // Fragment supports alt allele on BP1
  evidence.read1.isScanned = true;
  evidence.read2.isScanned = true;
  // Both Read-1 and Read-2 are anchored reads that means both read-1 and read-2
  // have confident mapping.
  // so power = spanningPairWeight = 0.2
  evidence.read1.setAnchored(true);
  evidence.read2.setAnchored(true);
  evidence.alt.bp1.isFragmentSupport = true;
  evidence.alt.bp1.fragLengthProb    = 0.4f;
  evidence.ref.bp1.fragLengthProb    = 0.6f;
  BOOST_REQUIRE(getRefAltFromFrag(
      spanningPairWeight,
      semiMappedPower,
      refChimeraProb,
      altChimeraProb,
      refSplitMapProb,
      altSplitMapProb,
      false,
      fragLabel,
      evidence,
      refLnLhoodSet,
      altLnLhoodSet,
      isRead1Evaluated,
      isRead2Evaluated));
  // Check for altLnLhood and refLnLhood
  BOOST_REQUIRE_CLOSE(refLnLhoodSet.fragPair, -0.2407945644533058, eps);
  BOOST_REQUIRE_CLOSE(altLnLhoodSet.fragPair, -0.14679383546496988, eps);

  // Fragment supports alt allele on BP1
  // Fragment supports ref allele on BP2
  // Here only read-1 has confident mapping. So it is semi mapped fragment
  // Also fragment-size probability of ref allele (0.6) is greater than alt allele(0.4)
  // so power = 0
  evidence.read2.setAnchored(false);
  evidence.ref.bp2.isFragmentSupport = true;
  BOOST_REQUIRE(getRefAltFromFrag(
      spanningPairWeight,
      semiMappedPower,
      refChimeraProb,
      altChimeraProb,
      refSplitMapProb,
      altSplitMapProb,
      false,
      fragLabel,
      evidence,
      refLnLhoodSet,
      altLnLhoodSet,
      isRead1Evaluated,
      isRead2Evaluated));
  // Check for altLnLhood and refLnLhood
  BOOST_REQUIRE_CLOSE(refLnLhoodSet.fragPair, -0.72238369335991737, eps);
  BOOST_REQUIRE_CLOSE(altLnLhoodSet.fragPair, -0.44038150639490964, eps);

  // Here only read-1 has confident mapping. So it is semi mapped fragment
  // Also fragment-size probability of alt allele (0.8) is greater than ref allele(0.6)
  // so power = spanningPairWeight * semiMappedPower = 0.4.
  evidence.alt.bp1.fragLengthProb = 0.8;
  BOOST_REQUIRE(getRefAltFromFrag(
      spanningPairWeight,
      semiMappedPower,
      refChimeraProb,
      altChimeraProb,
      refSplitMapProb,
      altSplitMapProb,
      false,
      fragLabel,
      evidence,
      refLnLhoodSet,
      altLnLhoodSet,
      isRead1Evaluated,
      isRead2Evaluated));
  // Check for altLnLhood and refLnLhood
  BOOST_REQUIRE_CLOSE(refLnLhoodSet.fragPair, -1.2039728222665289, eps);
  BOOST_REQUIRE_CLOSE(altLnLhoodSet.fragPair, -0.55015624191946355, eps);

  // Fragment does not support any allele on breakpoints
  evidence.read1.isScanned = false;
  evidence.read2.isScanned = false;
  evidence.read1.setAnchored(false);
  evidence.read2.setAnchored(false);
  evidence.alt.bp1.isFragmentSupport = false;
  // Only Read-1 supports split evidence on BP1. So according to case-1
  // api returns true.
  evidence.alt.bp1.read1.isSplitSupport   = true;
  evidence.alt.bp1.read1.isSplitEvaluated = true;
  evidence.alt.bp2.read1.isSplitEvaluated = true;
  evidence.ref.bp1.read1.isSplitEvaluated = true;
  evidence.ref.bp2.read1.isSplitEvaluated = true;
  BOOST_REQUIRE(getRefAltFromFrag(
      spanningPairWeight,
      semiMappedPower,
      refChimeraProb,
      altChimeraProb,
      refSplitMapProb,
      altSplitMapProb,
      false,
      fragLabel,
      evidence,
      refLnLhoodSet,
      altLnLhoodSet,
      isRead1Evaluated,
      isRead2Evaluated));
  BOOST_REQUIRE(isRead1Evaluated);
  BOOST_REQUIRE(!isRead2Evaluated);

  // Fragment does not support any allele on breakpoints.
  evidence.alt.bp1.read1.isSplitSupport = false;
  // Only Read-2 supports split evidence on BP1. So according to case-1
  // api returns true.
  evidence.alt.bp1.read2.isSplitSupport   = true;
  evidence.alt.bp1.read2.isSplitEvaluated = true;
  evidence.alt.bp2.read2.isSplitEvaluated = true;
  evidence.ref.bp1.read2.isSplitEvaluated = true;
  evidence.ref.bp2.read2.isSplitEvaluated = true;
  BOOST_REQUIRE(getRefAltFromFrag(
      spanningPairWeight,
      semiMappedPower,
      refChimeraProb,
      altChimeraProb,
      refSplitMapProb,
      altSplitMapProb,
      false,
      fragLabel,
      evidence,
      refLnLhoodSet,
      altLnLhoodSet,
      isRead1Evaluated,
      isRead2Evaluated));
  BOOST_REQUIRE(!isRead1Evaluated);
  BOOST_REQUIRE(isRead2Evaluated);
}

// Test the Diploid scoring model
// For reference and alternate alleles A = {r, x},
// the genotype states at each allele are restricted to G = {rr, rx, xx}.
// So if fragment-size likelihood of ref allele is R and fragment-size likelihood of alt allele
// is L, then score(S) of a genotype (gt = REF or HET or HOM) is
// S = log_sum(R + altLnCompFraction(gt), L + altLnFraction(gt)), where,
//     altLnFraction(gt) = log(Prior probability of expected alt allele)
//     altLnCompFraction(gt) = log(1-Prior probability of expected alt allele)
//     log_sum(x1,x2) = x1 + log(1+exp(x2-x1)) if x2>x1
//                    = x2 + log(1+exp(x1-x2)) if x1>x2
BOOST_AUTO_TEST_CASE(test_addDiploidLoglhood)
{
  static const double         eps(0.00000001);
  SVEvidence::evidenceTrack_t evidenceTrack;
  std::string                 fragLabel = "frag-1";
  SVFragmentEvidence          evidence;
  evidenceTrack[fragLabel] = evidence;
  std::array<double, DIPLOID_GT::SIZE> loglhood;
  for (unsigned gt(0); gt < DIPLOID_GT::SIZE; ++gt) {
    loglhood[gt] = 0;
  }
  // fragment was not evaluated for pair or split support for either allele
  addDiploidLoglhood(0.2, evidenceTrack, loglhood);
  for (unsigned gt(0); gt < DIPLOID_GT::SIZE; ++gt) {
    BOOST_REQUIRE(loglhood[gt] == 0);
  }
  // fragment was evaluated for alt allele
  evidence.read1.isScanned = true;
  evidence.read2.isScanned = true;
  evidence.read1.setAnchored(true);
  evidence.read2.setAnchored(true);
  evidence.alt.bp1.isFragmentSupport = true;
  evidenceTrack[fragLabel]           = evidence;
  addDiploidLoglhood(0.2, evidenceTrack, loglhood);
  // above mentioned formula is applied in loglhood
  BOOST_REQUIRE_CLOSE(loglhood[0], -1.3815510763831425, eps);  // REF
  BOOST_REQUIRE_CLOSE(loglhood[1], -1.3815510782877967, eps);  // HET
  BOOST_REQUIRE_CLOSE(loglhood[2], -1.3815510773843347, eps);  // HOM
}

// Test the following cases:
// 1. when quality of a diploid variant is less than minimum QUAL score (20). Add minAltFilterLabel.
// 2. If no sample passes all sample-specific filters, apply sample FT filter at the record level.
//    Add failedSampleFTLabel.
// 3. when MQ0 fraction of BP1 or BP2 is greater than control filtration based on MQ0 fraction (0.4). Add
//    maxMQ0FracLabel.
// 4. when half of the junctions do not support spanning pair confidently. Add noPairSupportLabel.
// 5. when max depth of BP1 or BP2 is greater than chromosome depth. Add maxDepthFilterLabel.
// 6. When genotype score less than threshold, add minGTFilterLabel.
// 7. When gentope is reference, add homRefLabel.
BOOST_AUTO_TEST_CASE(test_scoreDiploidSV)
{
  CallOptionsDiploid              diploidOpt;
  const bam_header_info           bamHeader(buildTestBamHeader());
  std::unique_ptr<SVLocusScanner> scanner(buildTestSVLocusScanner(bamHeader));
  CallOptionsDiploidDeriv         diploidDeriv(diploidOpt);
  // dummy depth file that means depth file is not present
  ChromDepthFilterUtil dFilter("", 2, bamHeader);

  // Diploid sample info
  SVScoreInfoDiploid diploidInfo;
  diploidInfo.setSampleCount(1);  // diploid sample count = 1
  SVFragmentEvidence fragmentEvidence;
  fragmentEvidence.read1.isScanned = true;
  fragmentEvidence.read2.isScanned = true;
  fragmentEvidence.read1.setAnchored(true);
  fragmentEvidence.read2.setAnchored(true);
  fragmentEvidence.alt.bp1.isFragmentSupport = true;
  SVEvidence::evidenceTrack_t evidenceTrack;
  std::string                 fragLabel = "frag-1";
  evidenceTrack[fragLabel]              = fragmentEvidence;
  SVEvidence evidence;
  evidence.samples.push_back(evidenceTrack);

  // Diploid score info
  SVScoreInfo  scoreInfo;
  SVSampleInfo diploidSampleInfo;
  diploidSampleInfo.alt.confidentSpanningPairCount           = 100;
  diploidSampleInfo.alt.confidentSemiMappedSpanningPairCount = 50;
  diploidSampleInfo.alt.confidentSplitReadCount              = 10;
  scoreInfo.samples.push_back(diploidSampleInfo);
  SVCandidate      candidate;
  JunctionCallInfo callInfo;
  callInfo.init(candidate, evidence, scoreInfo, 0.2);
  std::vector<JunctionCallInfo> junctionData = {callInfo};

  scoreDiploidSV(diploidOpt, scanner.operator*(), diploidDeriv, dFilter, junctionData, diploidInfo);
  // Designed the case-1 and case-2
  // Here phred score of the variant is 0 which is less than 20. Add a filter label minAltFilterLabel to this
  // record. Here we have only one diploid sample and this variant is reference genotype (based on maximum
  // value in loglnlhood in test_addDiploidLoglhood) for this sample. Then add failedSampleFTLabel to this
  // record.
  BOOST_REQUIRE_EQUAL(diploidInfo.filters.size(), 2);
  BOOST_REQUIRE_EQUAL(*(diploidInfo.filters.begin()), diploidOpt.minAltFilterLabel);
  BOOST_REQUIRE_EQUAL(*(std::next(diploidInfo.filters.begin(), 1)), diploidOpt.failedSampleFTLabel);

  junctionData.clear();
  diploidInfo.filters.clear();
  scoreInfo.bp1MQ0Frac = 0.2f;
  scoreInfo.bp2MQ0Frac = 0.5f;
  callInfo.init(candidate, evidence, scoreInfo, 0.2);
  junctionData.push_back(callInfo);
  scoreDiploidSV(diploidOpt, scanner.operator*(), diploidDeriv, dFilter, junctionData, diploidInfo);
  // Designed the case-1, case-2 and case-3
  // For case-3: If fraction of mapping quality 0 reads more than 0.4 at any of the breakpoints,
  // then add maxMQ0FracLabel label to this record.
  // Here fraction of mapping quality 0 reads at BP2 is 0.5 which is greater than 0.4.
  BOOST_REQUIRE_EQUAL(diploidInfo.filters.size(), 3);
  BOOST_REQUIRE_EQUAL(*(diploidInfo.filters.begin()), diploidOpt.maxMQ0FracLabel);
  BOOST_REQUIRE_EQUAL(*(std::next(diploidInfo.filters.begin(), 1)), diploidOpt.minAltFilterLabel);
  BOOST_REQUIRE_EQUAL(*(std::next(diploidInfo.filters.begin(), 2)), diploidOpt.failedSampleFTLabel);

  // Designed the case-1, case-2 and case-4
  candidate.bp1.interval = GenomeInterval(0, 100, 200);
  candidate.bp2.interval = GenomeInterval(1, 200, 300);
  junctionData.clear();
  diploidInfo.filters.clear();
  scoreInfo.samples[0].alt.confidentSpanningPairCount = 0;
  callInfo.init(candidate, evidence, scoreInfo, 0.2);
  junctionData.push_back(callInfo);
  // Let's say, a junction is called bad if there is no confident spanning pair support of alt allele for this
  // junction. If number of bad junction is more than half of the total junctions, then add
  // noPairSupportLabel to diploid info.
  // So here confidentSpanningPairCount of alt allele is 0 and total number of junction is 1.
  // So number of bad junction is also 1 which is greater than 1/2 = 0.5
  scoreDiploidSV(diploidOpt, scanner.operator*(), diploidDeriv, dFilter, junctionData, diploidInfo);
  BOOST_REQUIRE_EQUAL(diploidInfo.filters.size(), 3);
  BOOST_REQUIRE_EQUAL(*(diploidInfo.filters.begin()), diploidOpt.minAltFilterLabel);
  BOOST_REQUIRE_EQUAL(*(std::next(diploidInfo.filters.begin(), 1)), diploidOpt.noPairSupportLabel);
  BOOST_REQUIRE_EQUAL(*(std::next(diploidInfo.filters.begin(), 2)), diploidOpt.failedSampleFTLabel);

  // Designed the case-6 and case-7. MinGTScore = 50
  // Here genotype score is 48 which is less than 50.
  diploidOpt.minPassGTScore = 50;
  diploidInfo.filters.clear();
  scoreDiploidSV(diploidOpt, scanner.operator*(), diploidDeriv, dFilter, junctionData, diploidInfo);
  BOOST_REQUIRE_EQUAL(diploidInfo.samples[0].filters.size(), 2);
  BOOST_REQUIRE_EQUAL(*(diploidInfo.samples[0].filters.begin()), diploidOpt.homRefLabel);
  BOOST_REQUIRE_EQUAL(*(std::next(diploidInfo.samples[0].filters.begin(), 1)), diploidOpt.minGTFilterLabel);

  // Designed the case-5 where BP depth is more than chromosome depth
  const std::shared_ptr<TestChromosomeDepthFileMaker> depthFileMaker(new TestChromosomeDepthFileMaker());
  const std::string                                   depthFileName(depthFileMaker.operator*().getFilename());
  buildTestChromosomeDepthFile(depthFileName);
  ChromDepthFilterUtil depthFilterUtil2(depthFileName, 2, bamHeader);
  SVScoreInfo          scoreInfo2;
  scoreInfo2.bp2MaxDepth = 120;
  scoreInfo2.samples.resize(1);
  scoreInfo2.samples.push_back(diploidSampleInfo);
  junctionData.clear();
  diploidInfo.filters.clear();
  callInfo.init(candidate, evidence, scoreInfo2, 0.2);
  junctionData.push_back(callInfo);
  scoreDiploidSV(diploidOpt, scanner.operator*(), diploidDeriv, depthFilterUtil2, junctionData, diploidInfo);
  BOOST_REQUIRE_EQUAL(diploidInfo.filters.size(), 4);
  // Here max depth at BP2 = 120 which is greater than chrBar depth (16)
  BOOST_REQUIRE_EQUAL(*(diploidInfo.filters.begin()), diploidOpt.maxDepthFilterLabel);
  BOOST_REQUIRE_EQUAL(*(std::next(diploidInfo.filters.begin(), 1)), diploidOpt.minAltFilterLabel);
  BOOST_REQUIRE_EQUAL(*(std::next(diploidInfo.filters.begin(), 2)), diploidOpt.noPairSupportLabel);
  BOOST_REQUIRE_EQUAL(*(std::next(diploidInfo.filters.begin(), 3)), diploidOpt.failedSampleFTLabel);
}

// Test the following cases:
// 1. when MQ0 fraction of BP1 or BP2 is greater than control filtration based on MQ0 fraction (0.4)
// 2. when max depth of BP1 or BP2 is greater than chromosome depth
BOOST_AUTO_TEST_CASE(test_scoreTumorSV)
{
  const bam_header_info bamHeader(buildTestBamHeader());
  CallOptionsTumor      optionsTumor;
  // Dummy depth file that means depth file is not present
  ChromDepthFilterUtil depthFilterUtil1("", 0.2, bamHeader);
  SVCandidate          candidate;
  candidate.bp1.interval = GenomeInterval(0, 100, 200);
  candidate.bp2.interval = GenomeInterval(0, 400, 500);
  candidate.insertSeq    = "AGCTGACTGATCAGT";
  SVScoreInfoTumor scoreInfoTumor;
  SVScoreInfo      scoreInfo1;
  scoreInfo1.bp1MQ0Frac = 0.2f;
  scoreInfo1.bp2MQ0Frac = 0.3f;
  SVEvidence       evidence;
  JunctionCallInfo callInfo;
  callInfo.init(candidate, evidence, scoreInfo1, 0.2);
  std::vector<JunctionCallInfo> callInfos;
  callInfos.push_back(callInfo);
  scoreTumorSV(optionsTumor, depthFilterUtil1, callInfos, scoreInfoTumor);
  // All the cases are fine. No label will be added.
  BOOST_REQUIRE_EQUAL(scoreInfoTumor.filters.size(), 0);

  // For case-1: If fraction of mapping quality 0 reads more than 0.4 at any of the breakpoints,
  // then add maxMQ0FracLabel label to this record.
  // Here fraction of mapping quality 0 reads at BP2 is 0.5 which is greater than 0.4.
  callInfos.clear();
  scoreInfo1.bp1MQ0Frac = 0.2f;
  scoreInfo1.bp2MQ0Frac = 0.5f;
  callInfo.init(candidate, evidence, scoreInfo1, 0.2);
  callInfos.push_back(callInfo);
  scoreTumorSV(optionsTumor, depthFilterUtil1, callInfos, scoreInfoTumor);
  BOOST_REQUIRE_EQUAL(scoreInfoTumor.filters.size(), 1);
  BOOST_REQUIRE_EQUAL(*(scoreInfoTumor.filters.begin()), optionsTumor.maxMQ0FracLabel);

  // Designed the case-1 for BP1
  callInfos.clear();
  scoreInfoTumor.filters.clear();
  scoreInfo1.bp1MQ0Frac = 0.5f;
  scoreInfo1.bp2MQ0Frac = 0.2f;
  callInfo.init(candidate, evidence, scoreInfo1, 0.2);
  callInfos.push_back(callInfo);
  scoreTumorSV(optionsTumor, depthFilterUtil1, callInfos, scoreInfoTumor);
  BOOST_REQUIRE_EQUAL(scoreInfoTumor.filters.size(), 1);
  BOOST_REQUIRE_EQUAL(*(scoreInfoTumor.filters.begin()), optionsTumor.maxMQ0FracLabel);

  // Designed the case-1 for both BP1 and BP2
  callInfos.clear();
  scoreInfoTumor.filters.clear();
  scoreInfo1.bp1MQ0Frac = 0.5f;
  scoreInfo1.bp2MQ0Frac = 0.6f;
  callInfo.init(candidate, evidence, scoreInfo1, 0.2);
  callInfos.push_back(callInfo);
  scoreTumorSV(optionsTumor, depthFilterUtil1, callInfos, scoreInfoTumor);
  BOOST_REQUIRE_EQUAL(scoreInfoTumor.filters.size(), 1);
  BOOST_REQUIRE_EQUAL(*(scoreInfoTumor.filters.begin()), optionsTumor.maxMQ0FracLabel);

  // Designed the case-2 where BP depth is more than chromosome depth
  callInfos.clear();
  scoreInfoTumor.filters.clear();
  const std::shared_ptr<TestChromosomeDepthFileMaker> depthFileMaker(new TestChromosomeDepthFileMaker());
  const std::string                                   depthFileName(depthFileMaker.operator*().getFilename());
  buildTestChromosomeDepthFile(depthFileName);
  ChromDepthFilterUtil depthFilterUtil2(depthFileName, 2, bamHeader);
  SVScoreInfo          scoreInfo2;
  scoreInfo2.bp2MaxDepth = 120;
  callInfo.init(candidate, evidence, scoreInfo2, 0.2);
  callInfos.push_back(callInfo);
  // Here max depth at BP2 = 120 which is greater than chrBar depth (16)
  scoreTumorSV(optionsTumor, depthFilterUtil2, callInfos, scoreInfoTumor);
  BOOST_REQUIRE_EQUAL(scoreInfoTumor.filters.size(), 1);
  BOOST_REQUIRE_EQUAL(*(scoreInfoTumor.filters.begin()), optionsTumor.maxDepthFilterLabel);
}

// Following cases need to be tested:
// 1. When SV is imprecise. Add Imprecise label.
// 2. When difference between two breakpoints is less than Min length for passing fusions (100000).
//    Add a local label.
// 3. when split read count or confident spanning pair count is 0. Add a rnaFilterLabel.
BOOST_AUTO_TEST_CASE(test_scoreRNASV)
{
  SVSampleInfo   sampleInfo;
  SVCandidate    candidate;
  SVScoreInfoRna scoreInfoRna;
  scoreRNASV(sampleInfo, candidate, scoreInfoRna);
  // Designed the case-1 where candidate sv is imprecise.
  BOOST_REQUIRE_EQUAL(scoreInfoRna.filters.size(), 1);
  BOOST_REQUIRE_EQUAL(*(scoreInfoRna.filters.begin()), SVScoreInfoRna::impreciseLabel);

  // Designed the case-2 difference between two breakpoints is less than Min length for
  // passing fusions (100000)
  candidate.setPrecise();
  scoreInfoRna.filters.clear();
  // Here difference between two breakpoint is 425-150 = 275 < 100000
  candidate.bp1.interval                    = GenomeInterval(0, 100, 200);
  candidate.bp2.interval                    = GenomeInterval(0, 400, 450);
  sampleInfo.alt.splitReadCount             = 1;
  sampleInfo.alt.confidentSpanningPairCount = 1;
  scoreRNASV(sampleInfo, candidate, scoreInfoRna);
  BOOST_REQUIRE_EQUAL(scoreInfoRna.filters.size(), 1);
  BOOST_REQUIRE_EQUAL(*(scoreInfoRna.filters.begin()), SVScoreInfoRna::localLabel);

  // Designed the case-2 difference between two breakpoints is greater than Min length for
  // passing fusions (100000)
  scoreInfoRna.filters.clear();
  // Here difference between two breakpoint is 400025-150 = 399875 > 100000
  candidate.bp1.interval                    = GenomeInterval(0, 100, 200);
  candidate.bp2.interval                    = GenomeInterval(0, 400000, 400050);
  sampleInfo.alt.splitReadCount             = 1;
  sampleInfo.alt.confidentSpanningPairCount = 1;
  scoreRNASV(sampleInfo, candidate, scoreInfoRna);
  BOOST_REQUIRE_EQUAL(scoreInfoRna.filters.size(), 0);

  // Designed the case-3 when splitReadCount = 0 and confidentSpanningPairCount = 1
  scoreInfoRna.filters.clear();
  candidate.bp1.interval                    = GenomeInterval(0, 100, 200);
  candidate.bp2.interval                    = GenomeInterval(1, 400, 450);
  sampleInfo.alt.splitReadCount             = 0;
  sampleInfo.alt.confidentSpanningPairCount = 1;
  scoreRNASV(sampleInfo, candidate, scoreInfoRna);
  BOOST_REQUIRE_EQUAL(scoreInfoRna.filters.size(), 1);
  BOOST_REQUIRE_EQUAL(*(scoreInfoRna.filters.begin()), SVScoreInfoRna::rnaFilterLabel);

  // Designed the case-3 when splitReadCount = 1 and confidentSpanningPairCount = 0
  scoreInfoRna.filters.clear();
  candidate.bp1.interval                    = GenomeInterval(0, 100, 200);
  candidate.bp2.interval                    = GenomeInterval(1, 400, 450);
  sampleInfo.alt.splitReadCount             = 1;
  sampleInfo.alt.confidentSpanningPairCount = 0;
  scoreRNASV(sampleInfo, candidate, scoreInfoRna);
  BOOST_REQUIRE_EQUAL(scoreInfoRna.filters.size(), 1);
  BOOST_REQUIRE_EQUAL(*(scoreInfoRna.filters.begin()), SVScoreInfoRna::rnaFilterLabel);

  // Designed the case-3 when splitReadCount = 0 and confidentSpanningPairCount = 0
  scoreInfoRna.filters.clear();
  candidate.bp1.interval                    = GenomeInterval(0, 100, 200);
  candidate.bp2.interval                    = GenomeInterval(0, 400, 450);
  sampleInfo.alt.splitReadCount             = 0;
  sampleInfo.alt.confidentSpanningPairCount = 0;
  scoreRNASV(sampleInfo, candidate, scoreInfoRna);
  BOOST_REQUIRE_EQUAL(scoreInfoRna.filters.size(), 2);
  BOOST_REQUIRE_EQUAL(*(scoreInfoRna.filters.begin()), SVScoreInfoRna::localLabel);
  BOOST_REQUIRE_EQUAL(*(std::next(scoreInfoRna.filters.begin(), 1)), SVScoreInfoRna::rnaFilterLabel);
}

// If tier2 support
//     spanning pair count = confidentSemiMappedSpanningPairCount * spanning pair weight
// else
//     spanning pair count = confidentSpanningPairCount * spanning pair weight
// where spanning pair weight is 0.2 and isTier2/isPermissive is the 3rd arg of
// function getSpanningPairCount().
BOOST_AUTO_TEST_CASE(test_getSpanningPairCount)
{
  SVSampleAlleleInfo alleleInfo;
  alleleInfo.confidentSpanningPairCount           = 100;
  alleleInfo.confidentSemiMappedSpanningPairCount = 50;
  BOOST_REQUIRE_EQUAL(getSpanningPairCount(alleleInfo, 0.2, false), 20);
  // for tier2 support
  BOOST_REQUIRE_EQUAL(getSpanningPairCount(alleleInfo, 0.2, true), 10);
}

// If tier2 support
//     spanning pair count = confidentSplitReadCount + confidentSemiMappedSpanningPairCount * spanning pair
//     weight
// else
//     spanning pair count = confidentSplitReadCount + confidentSpanningPairCount * spanning pair weight
// where spanning pair weight is 0.2
BOOST_AUTO_TEST_CASE(test_getSupportCount)
{
  SVSampleAlleleInfo alleleInfo;
  alleleInfo.confidentSpanningPairCount           = 100;
  alleleInfo.confidentSemiMappedSpanningPairCount = 50;
  alleleInfo.confidentSplitReadCount              = 10;
  BOOST_REQUIRE_EQUAL(getSupportCount(alleleInfo, 0.2, false), 30);
  BOOST_REQUIRE_EQUAL(getSupportCount(alleleInfo, 0.2, true), 20);
}

// Test the fraction of allele support count to total support
// counts (alt support + ref support)
BOOST_AUTO_TEST_CASE(test_estimateSomaticMutationFreq)
{
  SVCandidate candidate;
  SVEvidence  evidence;

  // All counts are zero so fraction should be also 0.
  std::vector<JunctionCallInfo> callInfos;
  callInfos.resize(1);
  SVScoreInfo  scoreInfo1;
  SVSampleInfo sampleInfo;
  scoreInfo1.samples.resize(1);
  scoreInfo1.samples.push_back(sampleInfo);
  callInfos[0].init(candidate, evidence, scoreInfo1, 0.2);
  BOOST_REQUIRE_EQUAL(estimateSomaticMutationFreq(0, callInfos, false), 0);

  // alt support count > 0 and ref support count is 0. So fraction of
  // alt support = 1
  sampleInfo.alt.confidentSpanningPairCount           = 100;
  sampleInfo.alt.confidentSemiMappedSpanningPairCount = 50;
  sampleInfo.alt.confidentSplitReadCount              = 10;
  SVScoreInfo scoreInfo2;
  scoreInfo2.samples.push_back(sampleInfo);
  callInfos[0].init(candidate, evidence, scoreInfo2, 0.2);
  BOOST_REQUIRE_EQUAL(estimateSomaticMutationFreq(0, callInfos, false), 1);

  // alt support count = 0 and ref support count > 0. So fraction of
  // alt support = 0
  sampleInfo.alt.confidentSpanningPairCount           = 0;
  sampleInfo.alt.confidentSemiMappedSpanningPairCount = 0;
  sampleInfo.alt.confidentSplitReadCount              = 0;
  sampleInfo.ref.confidentSpanningPairCount           = 100;
  sampleInfo.ref.confidentSemiMappedSpanningPairCount = 50;
  sampleInfo.ref.confidentSplitReadCount              = 10;
  SVScoreInfo scoreInfo3;
  scoreInfo3.samples.push_back(sampleInfo);
  callInfos[0].init(candidate, evidence, scoreInfo3, 0.2);
  BOOST_REQUIRE_EQUAL(estimateSomaticMutationFreq(0, callInfos, false), 0);

  // alt support count > 0 and ref support count is >. So fraction of
  // alt support = alt support / (alt support + ref support)
  sampleInfo.alt.confidentSpanningPairCount           = 100;
  sampleInfo.alt.confidentSemiMappedSpanningPairCount = 50;
  sampleInfo.alt.confidentSplitReadCount              = 10;
  sampleInfo.ref.confidentSpanningPairCount           = 100;
  sampleInfo.ref.confidentSemiMappedSpanningPairCount = 50;
  sampleInfo.ref.confidentSplitReadCount              = 10;
  SVScoreInfo scoreInfo4;
  scoreInfo4.samples.push_back(sampleInfo);
  callInfos[0].init(candidate, evidence, scoreInfo4, 0.2);
  BOOST_REQUIRE_EQUAL(estimateSomaticMutationFreq(0, callInfos, false), 0.5);
}

// This is applicable for somatic caller.
// Test the fraction of somatic alt support count (normal + tumor) to total support
// counts (alt support + ref support of normal and tumor)
BOOST_AUTO_TEST_CASE(test_estimateNoiseMutationFreq)
{
  SVCandidate candidate;
  SVEvidence  evidence;
  // All counts are zero for both normal and tumor
  std::vector<JunctionCallInfo> callInfos;
  callInfos.resize(1);
  SVScoreInfo  scoreInfo1;
  SVSampleInfo normalSampleInfo;
  SVSampleInfo tumorSampleInfo;
  scoreInfo1.samples.resize(2);
  scoreInfo1.samples.push_back(normalSampleInfo);
  scoreInfo1.samples.push_back(tumorSampleInfo);
  callInfos[0].init(candidate, evidence, scoreInfo1, 0.2);
  BOOST_REQUIRE_EQUAL(estimateNoiseMutationFreq(0, 1, callInfos, false), 0);

  // alt support count (tumor + normal) > 0 and ref support count is 0. So fraction of somatic
  // alt support = 1
  normalSampleInfo.alt.confidentSpanningPairCount           = 100;
  normalSampleInfo.alt.confidentSemiMappedSpanningPairCount = 50;
  normalSampleInfo.alt.confidentSplitReadCount              = 10;
  tumorSampleInfo.alt.confidentSpanningPairCount            = 100;
  tumorSampleInfo.alt.confidentSemiMappedSpanningPairCount  = 50;
  tumorSampleInfo.alt.confidentSplitReadCount               = 10;
  SVScoreInfo scoreInfo2;
  scoreInfo2.samples.push_back(normalSampleInfo);
  scoreInfo2.samples.push_back(tumorSampleInfo);
  callInfos[0].init(candidate, evidence, scoreInfo2, 0.2);
  BOOST_REQUIRE_EQUAL(estimateNoiseMutationFreq(0, 1, callInfos, false), 1);

  // alt support count (tumor + normal) = 0 and ref support count
  // (tumor + normal) > 0. So fraction of somatic alt support = 0
  normalSampleInfo.alt.confidentSpanningPairCount           = 0;
  normalSampleInfo.alt.confidentSemiMappedSpanningPairCount = 0;
  normalSampleInfo.alt.confidentSplitReadCount              = 0;
  normalSampleInfo.ref.confidentSpanningPairCount           = 100;
  normalSampleInfo.ref.confidentSemiMappedSpanningPairCount = 50;
  normalSampleInfo.ref.confidentSplitReadCount              = 10;
  tumorSampleInfo.alt.confidentSpanningPairCount            = 0;
  tumorSampleInfo.alt.confidentSemiMappedSpanningPairCount  = 0;
  tumorSampleInfo.alt.confidentSplitReadCount               = 0;
  tumorSampleInfo.ref.confidentSpanningPairCount            = 100;
  tumorSampleInfo.ref.confidentSemiMappedSpanningPairCount  = 50;
  tumorSampleInfo.ref.confidentSplitReadCount               = 10;
  SVScoreInfo scoreInfo3;
  scoreInfo3.samples.push_back(normalSampleInfo);
  scoreInfo3.samples.push_back(tumorSampleInfo);
  callInfos[0].init(candidate, evidence, scoreInfo3, 0.2);
  BOOST_REQUIRE_EQUAL(estimateNoiseMutationFreq(0, 1, callInfos, false), 0);

  // alt support count (tumor + normal) > 0 and ref support count (tumor + normal) > 0.
  // So fraction of somatic alt support = alt support / (alt support + ref support)
  normalSampleInfo.alt.confidentSpanningPairCount           = 100;
  normalSampleInfo.alt.confidentSemiMappedSpanningPairCount = 50;
  normalSampleInfo.alt.confidentSplitReadCount              = 10;
  normalSampleInfo.ref.confidentSpanningPairCount           = 100;
  normalSampleInfo.ref.confidentSemiMappedSpanningPairCount = 50;
  normalSampleInfo.ref.confidentSplitReadCount              = 10;
  tumorSampleInfo.alt.confidentSpanningPairCount            = 100;
  tumorSampleInfo.alt.confidentSemiMappedSpanningPairCount  = 50;
  tumorSampleInfo.alt.confidentSplitReadCount               = 10;
  tumorSampleInfo.ref.confidentSpanningPairCount            = 100;
  tumorSampleInfo.ref.confidentSemiMappedSpanningPairCount  = 50;
  tumorSampleInfo.ref.confidentSplitReadCount               = 10;
  SVScoreInfo scoreInfo4;
  scoreInfo4.samples.push_back(normalSampleInfo);
  scoreInfo4.samples.push_back(tumorSampleInfo);
  callInfos[0].init(candidate, evidence, scoreInfo4, 0.2);
  BOOST_REQUIRE_EQUAL(estimateNoiseMutationFreq(0, 1, callInfos, false), 0.5);
}

// The simplified form a somatic genotype state space
// is defined consisting of non-somatic germline variant states {rr, rx, xx}, a noise
// state n representing spurious observations at the same allele frequency in the
// tumor and normal sample and the somatic state s.
// So if fragment-size likelihood of ref allele is R and fragment-size likelihood of alt allele
// is L, then score(S) of a genotype (gt = REF or HET or HOM or n or s) is
// S = log_sum(R + altLnCompFraction(gt), L + altLnFraction(gt)), where,
//     altLnFraction(gt) = log(Prior probability of expected alt allele)
//     altLnCompFraction(gt) = log(1-Prior probability of expected alt allele)
//     log_sum(x1,x2) = x1 + log(1+exp(x2-x1)) if x2>x1
//                    = x2 + log(1+exp(x1-x2)) if x1>x2
BOOST_AUTO_TEST_CASE(test_computeSomaticSampleLoghood)
{
  static const double         eps(0.00000001);
  float                       spanningPairWeight(0.2);
  ProbSet                     refChimeraProb(0.2);
  ProbSet                     altChimeraProb(0.3);
  ProbSet                     refSplitMapProb(0.1);
  ProbSet                     altSplitMapProb(0.5);
  SVEvidence::evidenceTrack_t evidenceTrack;
  std::string                 fragLabel = "frag-1";
  SVFragmentEvidence          fragmentEvidence1;
  fragmentEvidence1.alt.bp1.isFragmentSupport    = true;
  fragmentEvidence1.alt.bp1.fragLengthProb       = 0.4;
  fragmentEvidence1.ref.bp1.read1.isSplitSupport = true;
  fragmentEvidence1.ref.bp1.read1.splitLnLhood   = 0.2;
  evidenceTrack[fragLabel]                       = fragmentEvidence1;
  std::array<double, SOMATIC_GT::SIZE> loglhood;
  for (unsigned gt(0); gt < SOMATIC_GT::SIZE; ++gt) {
    loglhood[gt] = 0;
  }
  // fragment was not evaluated for pair or split support for either allele
  computeSomaticSampleLoghood(
      spanningPairWeight,
      evidenceTrack,
      0.4,
      0.25,
      false,
      false,
      refChimeraProb,
      altChimeraProb,
      refSplitMapProb,
      altSplitMapProb,
      loglhood);
  for (unsigned gt(0); gt < SOMATIC_GT::SIZE; ++gt) {
    BOOST_REQUIRE_EQUAL(loglhood[gt], 0);
  }

  // fragment was evaluated for alt allele
  fragmentEvidence1.read1.isScanned = true;
  fragmentEvidence1.read2.isScanned = true;
  fragmentEvidence1.read1.setAnchored(true);
  fragmentEvidence1.read2.setAnchored(true);
  fragmentEvidence1.alt.bp1.isFragmentSupport = true;
  evidenceTrack[fragLabel]                    = fragmentEvidence1;
  computeSomaticSampleLoghood(
      spanningPairWeight,
      evidenceTrack,
      0.4,
      0.25,
      false,
      false,
      refChimeraProb,
      altChimeraProb,
      refSplitMapProb,
      altSplitMapProb,
      loglhood);
  // above mentioned formula is applied in loglhood
  BOOST_REQUIRE_CLOSE(loglhood[0], -0.2407945644533058, eps);   // REF
  BOOST_REQUIRE_CLOSE(loglhood[1], -0.19269008924115461, eps);  // HET
  BOOST_REQUIRE_CLOSE(loglhood[2], -0.14679383546496988, eps);  // HOM
  BOOST_REQUIRE_CLOSE(loglhood[3], -0.20212764047652254, eps);  // SOM
  BOOST_REQUIRE_CLOSE(loglhood[4], -0.21645309966549675, eps);  // NOISE
}

// Following cases need to be tested
// 1. When somatic variant score is less than minimum somatic quality(30). Add minScoreLabel.
// 2. when max depth of BP1 or BP2 is greater than chromosome depth. Add maxDepth Label.
// 3. when MQ0 fraction of BP1 or BP2 is greater than control filtration based on MQ0 fraction (0.4). Add
//    maxMQ0FracLabel
BOOST_AUTO_TEST_CASE(test_scoreSomaticSV)
{
  // somatic caller options
  CallOptionsSomatic      somaticOpt;
  CallOptionsSomaticDeriv somaticDopt(somaticOpt);
  // created chromsome depth file
  const std::shared_ptr<TestChromosomeDepthFileMaker> depthFileMaker(new TestChromosomeDepthFileMaker());
  const std::string                                   depthFileName(depthFileMaker.operator*().getFilename());
  buildTestChromosomeDepthFile(depthFileName);
  const bam_header_info bamHeader(buildTestBamHeader());
  ChromDepthFilterUtil  depthFilterUtil(depthFileName, 2, bamHeader);

  std::vector<JunctionCallInfo> junctionData;
  junctionData.resize(1);  // adding one junction
  // Adding few counts to tumor and normal sample
  SVSampleInfo normalSampleInfo;
  SVSampleInfo tumorSampleInfo;
  normalSampleInfo.alt.confidentSpanningPairCount           = 100;
  normalSampleInfo.alt.confidentSemiMappedSpanningPairCount = 50;
  normalSampleInfo.alt.confidentSplitReadCount              = 10;
  normalSampleInfo.ref.confidentSpanningPairCount           = 100;
  normalSampleInfo.ref.confidentSemiMappedSpanningPairCount = 50;
  normalSampleInfo.ref.confidentSplitReadCount              = 10;
  tumorSampleInfo.alt.confidentSpanningPairCount            = 100;
  tumorSampleInfo.alt.confidentSemiMappedSpanningPairCount  = 50;
  tumorSampleInfo.alt.confidentSplitReadCount               = 10;
  tumorSampleInfo.ref.confidentSpanningPairCount            = 100;
  tumorSampleInfo.ref.confidentSemiMappedSpanningPairCount  = 50;
  tumorSampleInfo.ref.confidentSplitReadCount               = 10;
  // somatic SV score
  SVScoreInfo scoreInfo;
  scoreInfo.samples.push_back(normalSampleInfo);
  scoreInfo.samples.push_back(tumorSampleInfo);
  SVEvidence evidence;
  evidence.samples.resize(2);
  SVCandidate candidate;
  junctionData[0].init(candidate, evidence, scoreInfo, 0.2);
  SVScoreInfoSomatic somaticInfo;
  // Designed the case-1 where somatic score is 0 which is less than 30,
  // then add minSomaticScoreLabel label to this record.
  scoreSomaticSV(2, 1, somaticOpt, somaticDopt, depthFilterUtil, junctionData, somaticInfo);
  BOOST_REQUIRE_EQUAL(somaticInfo.filters.size(), 1);
  BOOST_REQUIRE_EQUAL(*(somaticInfo.filters.begin()), somaticOpt.minSomaticScoreLabel);

  // Designed the case-1 and case-2
  // For case-2, max depth at BP2 = 120 which is greater than chrBar depth (16).
  // Add maxDepthFilterLabel to this record.
  scoreInfo.bp2MaxDepth = 120;
  somaticInfo.filters.clear();
  junctionData[0].init(candidate, evidence, scoreInfo, 0.2);
  scoreSomaticSV(2, 1, somaticOpt, somaticDopt, depthFilterUtil, junctionData, somaticInfo);
  BOOST_REQUIRE_EQUAL(somaticInfo.filters.size(), 2);
  BOOST_REQUIRE_EQUAL(*(somaticInfo.filters.begin()), somaticOpt.maxDepthFilterLabel);
  BOOST_REQUIRE_EQUAL(*(std::next(somaticInfo.filters.begin(), 1)), somaticOpt.minSomaticScoreLabel);

  // Designed the case-1 and case-3
  // For case-3, If fraction of mapping quality 0 reads more than 0.4 at any of the breakpoints,
  // then add maxMQ0FracLabel label to this record.
  // Here fraction of mapping quality 0 reads at BP2 is 0.6 which is greater than 0.4.
  scoreInfo.bp2MaxDepth = 20;
  scoreInfo.bp2MQ0Frac  = 0.6;
  somaticInfo.filters.clear();
  somaticInfo.somaticScore = 40;
  junctionData[0].init(candidate, evidence, scoreInfo, 0.2);
  scoreSomaticSV(2, 1, somaticOpt, somaticDopt, depthFilterUtil, junctionData, somaticInfo);
  BOOST_REQUIRE_EQUAL(somaticInfo.filters.size(), 2);
  BOOST_REQUIRE_EQUAL(*(somaticInfo.filters.begin()), somaticOpt.maxMQ0FracLabel);
  BOOST_REQUIRE_EQUAL(*(std::next(somaticInfo.filters.begin(), 1)), somaticOpt.minSomaticScoreLabel);
}

// Test the fraction of zero mapping quality reads
BOOST_AUTO_TEST_CASE(test_getBreakendMaxMappedDepthAndMQ0)
{
  static const double             eps(0.00000001);
  const bam_header_info           bamHeader(buildTestBamHeader());
  std::unique_ptr<SVLocusScanner> scanner(buildTestSVLocusScanner(bamHeader));
  GSCOptions                      options;
  options.alignFileOpt.alignmentFilenames = {bamFilename};
  options.alignFileOpt.isAlignmentTumor   = {false};
  SVScorer     scorer(options, scanner.operator*(), bamHeader);
  TestSVScorer fSVScorer;
  SVBreakend   breakend;
  breakend.interval = GenomeInterval(0, 109, 110);
  breakend.state    = SVBreakendState::RIGHT_OPEN;
  float    mq0Frac(0);
  unsigned maxDepth(1000);
  // No depth cutoff in this test cases
  // Total 18 reads have mapping quality 0 out of 20 reads. So fraction is 18/20 = 0.9
  fSVScorer.getBreakendMaxMappedDepthAndMQ0(scorer, false, false, 100, breakend, maxDepth, mq0Frac);
  BOOST_REQUIRE_CLOSE(mq0Frac, 0.9f, eps);

  // depth cutoff is 12. If depth of a location more than 12, it will return the fraction.
  // So here total 11 reads have mapping quality 0 out of 13 reads. So fraction is 11/13 = ~0.85.
  fSVScorer.getBreakendMaxMappedDepthAndMQ0(scorer, false, true, 12, breakend, maxDepth, mq0Frac);
  BOOST_REQUIRE_CLOSE(mq0Frac, 0.846153855f, eps);
}

// Test Whether a specific scoring model is performed based on arguments
BOOST_AUTO_TEST_CASE(test_computeAllScoreModel)
{
  bool                            isSomatic(true);
  bool                            isTumorOnly(true);
  const bam_header_info           bamHeader(buildTestBamHeader());
  std::unique_ptr<SVLocusScanner> scanner(buildTestSVLocusScanner(bamHeader));
  GSCOptions                      options;
  options.alignFileOpt.alignmentFilenames = {bamFilename};
  SVScorer         scorer(options, scanner.operator*(), bamHeader);
  TestSVScorer     fSVScorer;
  SVModelScoreInfo modelScoreInfo;
  modelScoreInfo.setSampleCount(1, 1);

  SVCandidate candidate;
  candidate.bp1.interval = GenomeInterval(0, 100, 200);
  candidate.bp2.interval = GenomeInterval(0, 400, 500);
  candidate.insertSeq    = "AGCTGACTGATCAGT";
  SVScoreInfo scoreInfo;
  scoreInfo.bp1MQ0Frac = 0.2f;
  scoreInfo.bp2MQ0Frac = 0.5f;
  SVEvidence       evidence;
  JunctionCallInfo callInfo;
  callInfo.init(candidate, evidence, scoreInfo, 0.2);
  std::vector<JunctionCallInfo> callInfos;
  callInfos.push_back(callInfo);
  // Tumor scoring model
  fSVScorer.computeAllScoreModels(scorer, !isSomatic, isTumorOnly, callInfos, modelScoreInfo);
  BOOST_REQUIRE_EQUAL(modelScoreInfo.tumor.filters.size(), 1);
  callInfos.clear();
  candidate.setPrecise();
  candidate.bp1.interval                                        = GenomeInterval(0, 100, 200);
  candidate.bp2.interval                                        = GenomeInterval(0, 400, 450);
  modelScoreInfo.base.samples[0].alt.splitReadCount             = 1;
  modelScoreInfo.base.samples[0].alt.confidentSpanningPairCount = 1;
  SVFragmentEvidence fragmentEvidence;
  fragmentEvidence.read1.isScanned = true;
  fragmentEvidence.read2.isScanned = true;
  fragmentEvidence.read1.setAnchored(true);
  fragmentEvidence.read2.setAnchored(true);
  fragmentEvidence.alt.bp1.isFragmentSupport = true;
  SVEvidence::evidenceTrack_t evidenceTrack;
  std::string                 fragLabel = "frag-1";
  evidenceTrack[fragLabel]              = fragmentEvidence;
  evidence.samples.push_back(evidenceTrack);
  scoreInfo.samples.resize(1);
  callInfo.init(candidate, evidence, scoreInfo, 0.2);
  callInfos.push_back(callInfo);

  // Only Diploid scoring model
  fSVScorer.computeAllScoreModels(scorer, !isSomatic, !isTumorOnly, callInfos, modelScoreInfo);
  BOOST_REQUIRE_EQUAL(modelScoreInfo.diploid.filters.size(), 4);

  // somatic and diploid scoring model
  CallOptionsSomatic                                  somaticOpt;
  CallOptionsSomaticDeriv                             somaticDopt(somaticOpt);
  const std::shared_ptr<TestChromosomeDepthFileMaker> depthFileMaker(new TestChromosomeDepthFileMaker());
  const std::string                                   depthFileName(depthFileMaker.operator*().getFilename());
  buildTestChromosomeDepthFile(depthFileName);
  ChromDepthFilterUtil          depthFilterUtil(depthFileName, 2, bamHeader);
  std::vector<JunctionCallInfo> junctionData;
  junctionData.resize(1);
  SVSampleInfo normalSampleInfo;
  SVSampleInfo tumorSampleInfo;
  normalSampleInfo.alt.confidentSpanningPairCount           = 100;
  normalSampleInfo.alt.confidentSemiMappedSpanningPairCount = 50;
  normalSampleInfo.alt.confidentSplitReadCount              = 10;
  normalSampleInfo.ref.confidentSpanningPairCount           = 100;
  normalSampleInfo.ref.confidentSemiMappedSpanningPairCount = 50;
  normalSampleInfo.ref.confidentSplitReadCount              = 10;
  tumorSampleInfo.alt.confidentSpanningPairCount            = 100;
  tumorSampleInfo.alt.confidentSemiMappedSpanningPairCount  = 50;
  tumorSampleInfo.alt.confidentSplitReadCount               = 10;
  tumorSampleInfo.ref.confidentSpanningPairCount            = 100;
  tumorSampleInfo.ref.confidentSemiMappedSpanningPairCount  = 50;
  tumorSampleInfo.ref.confidentSplitReadCount               = 10;
  SVScoreInfo scoreInfo2;
  scoreInfo2.samples.push_back(normalSampleInfo);
  scoreInfo2.samples.push_back(tumorSampleInfo);
  SVEvidence evidence2;
  evidence2.samples.resize(2);
  SVCandidate candidate2;
  junctionData[0].init(candidate2, evidence2, scoreInfo2, 0.2);
  SVScoreInfoSomatic somaticInfo;
  somaticInfo.somaticScore = 20;
  modelScoreInfo.setSampleCount(2, 1);
  fSVScorer.computeAllScoreModels(scorer, isSomatic, !isTumorOnly, junctionData, modelScoreInfo);
  BOOST_REQUIRE_EQUAL(modelScoreInfo.somatic.filters.size(), 1);

  // RNA scoring model
  modelScoreInfo.setSampleCount(1, 1);
  GSCOptions options2;
  options2.alignFileOpt.alignmentFilenames = {bamFilename};
  options2.isRNA                           = true;
  SVScorer     scorer2(options2, scanner.operator*(), bamHeader);
  TestSVScorer fSVScorer2;
  fSVScorer2.computeAllScoreModels(scorer2, !isSomatic, !isTumorOnly, callInfos, modelScoreInfo);
  BOOST_REQUIRE_EQUAL(modelScoreInfo.rna.filters.size(), 2);
}

// Few minor logic checking have been done in this test case
// 1. When no valid junction is present. It will throw an exception
// 2. When only one unfiltered junction is present
// 3. When more than one unfiltered junction is present
BOOST_AUTO_TEST_CASE(test_ScoreSV)
{
  // Designed case-1 where there is no sv junction
  // It will throw an exception
  SVCandidateSetData                   svData;
  std::vector<SVCandidateAssemblyData> mjAssemblyData;
  SVMultiJunctionCandidate             mjSV;
  std::vector<SVId>                    mjSVId;
  std::vector<bool>                    isJunctionFiltered;
  std::vector<SVModelScoreInfo>        mjModelScoreInfo;
  SVModelScoreInfo                     mjJointModelScoreInfo;
  bool                                 isMJEvent;
  SVEvidenceWriterData                 svEvidenceWriterData(1);
  const bam_header_info                bamHeader(buildTestBamHeader());
  std::unique_ptr<SVLocusScanner>      scanner(buildTestSVLocusScanner(bamHeader));
  GSCOptions                           options;
  options.alignFileOpt.alignmentFilenames = {bamFilename};
  options.alignFileOpt.isAlignmentTumor   = {false};
  SVScorer scorer(options, scanner.operator*(), bamHeader);
  BOOST_CHECK_THROW(
      scorer.scoreSV(
          svData,
          mjAssemblyData,
          mjSV,
          mjSVId,
          isJunctionFiltered,
          false,
          true,
          mjModelScoreInfo,
          mjJointModelScoreInfo,
          isMJEvent,
          svEvidenceWriterData),
      illumina::common::GeneralException);

  // Designed case-2. Here number of unfiltered multi-junction event is 1.
  // So isMJEvent should be false here.
  SVCandidate candidate1;
  candidate1.bp1.interval = GenomeInterval(0, 100, 200);
  candidate1.bp2.interval = GenomeInterval(0, 300, 350);
  candidate1.bp1.state    = SVBreakendState::RIGHT_OPEN;
  candidate1.bp2.state    = SVBreakendState::LEFT_OPEN;
  candidate1.insertSeq    = "AGCTGACTGATCAGT";
  candidate1.setPrecise();
  candidate1.assemblyAlignIndex = 0;
  mjSV.junction.push_back(candidate1);
  SVId id1;
  mjSVId.push_back(id1);
  isJunctionFiltered.push_back(false);
  SVModelScoreInfo modelScoreInfo;
  modelScoreInfo.setSampleCount(1, 1);
  mjModelScoreInfo.push_back(modelScoreInfo);
  mjJointModelScoreInfo.setSampleCount(1, 1);
  TestSVScorer fSVScorer;
  fSVScorer.setSampleCount(scorer, 1, 1);

  std::string queryseq1 = "AGCTGACTGATCGATTTTTTACGTAGAGGAGCTTTGACGTATGAGCCTGATATGAGCCTG";
  std::string queryseq2 = "TGACGTATGAGCCTGATATGAGCCT";
  bam_record  bamRecord1;
  buildTestBamRecord(bamRecord1, 0, 200, 0, 300, 125, 15, "60M", queryseq1);
  bamRecord1.set_qname("Read-1");
  bamRecord1.toggle_is_first();

  bam_record bamRecord2;
  buildTestBamRecord(bamRecord2, 0, 300, 0, 200, 125, 15, "25M", queryseq2);
  bamRecord2.toggle_is_second();
  bamRecord2.toggle_is_fwd_strand();
  bamRecord2.toggle_is_mate_fwd_strand();
  bamRecord2.set_qname("Read-1");
  SVCandidateSetSequenceFragmentSampleGroup& group = svData.getDataGroup(0);
  group.add(bamHeader, bamRecord1, false, true, true);
  group.add(bamHeader, bamRecord2, false, true, true);
  SVSequenceFragmentAssociation association(0, SVEvidenceType::PAIR);
  group.begin().                operator*().svLink.push_back(association);

  SVCandidateAssemblyData candidateAssemblyData1;
  candidateAssemblyData1.bestAlignmentIndex = 0;
  Alignment alignment1;
  alignment1.beginPos = 406;
  std::string testCigar1("94=");
  cigar_to_apath(testCigar1.c_str(), alignment1.apath);
  Alignment alignment2;
  alignment2.beginPos = 510;
  std::string testCigar2("110=75I1D2=");
  cigar_to_apath(testCigar2.c_str(), alignment2.apath);
  SVCandidateAssemblyData::JumpAlignmentResultType jumpAlignmentResultType1;
  jumpAlignmentResultType1.align1    = alignment1;
  jumpAlignmentResultType1.align2    = alignment2;
  jumpAlignmentResultType1.jumpRange = 2;
  candidateAssemblyData1.spanningAlignments.push_back(jumpAlignmentResultType1);
  candidateAssemblyData1.isSpanning = true;
  candidateAssemblyData1.extendedContigs.push_back(
      "GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGT"
      "ATTTTCGTCTGGGGGGTGTGCACGCGATAGCATTGCGAGACGCTGGA");
  mjAssemblyData.push_back(candidateAssemblyData1);
  scorer.scoreSV(
      svData,
      mjAssemblyData,
      mjSV,
      mjSVId,
      isJunctionFiltered,
      false,
      true,
      mjModelScoreInfo,
      mjJointModelScoreInfo,
      isMJEvent,
      svEvidenceWriterData);
  BOOST_REQUIRE(!isMJEvent);

  // Designed case-3. Here number of unfiltered multi-junction event is 2.
  // So isMJEvent should be true here.
  SVCandidate candidate2;
  candidate2.bp1.interval = GenomeInterval(0, 400, 410);
  candidate2.bp2.interval = GenomeInterval(0, 505, 520);
  candidate2.bp1.state    = SVBreakendState::RIGHT_OPEN;
  candidate2.bp2.state    = SVBreakendState::LEFT_OPEN;
  candidate2.insertSeq    = "AGCTGACTGATCAGT";
  candidate2.setPrecise();
  candidate2.assemblyAlignIndex = 0;
  mjSV.junction.push_back(candidate2);
  SVId id2;
  mjSVId.push_back(id2);
  isJunctionFiltered.push_back(false);
  modelScoreInfo.setSampleCount(1, 1);
  mjModelScoreInfo.push_back(modelScoreInfo);
  mjJointModelScoreInfo.setSampleCount(1, 1);
  SVCandidateAssemblyData candidateAssemblyData2;
  candidateAssemblyData2.bestAlignmentIndex = 0;
  Alignment alignment3;
  alignment3.beginPos = 406;
  std::string testCigar3("54=");
  cigar_to_apath(testCigar3.c_str(), alignment3.apath);
  Alignment alignment4;
  alignment4.beginPos = 510;
  std::string testCigar4("110=75I1D2=");
  cigar_to_apath(testCigar4.c_str(), alignment4.apath);
  SVCandidateAssemblyData::JumpAlignmentResultType jumpAlignmentResultType2;
  jumpAlignmentResultType2.align1    = alignment3;
  jumpAlignmentResultType2.align2    = alignment4;
  jumpAlignmentResultType2.jumpRange = 2;
  candidateAssemblyData2.spanningAlignments.push_back(jumpAlignmentResultType2);
  candidateAssemblyData2.isSpanning = true;
  candidateAssemblyData2.extendedContigs.push_back(
      "GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGT"
      "ATTTTCGTCTGGGGGGTGTGCACGCGATAGCATTGCGAGACGCTGGA");
  mjAssemblyData.push_back(candidateAssemblyData2);
  //SVScorer scorer2(options, buildSVLocusScannerForSomatic(bamHeader).operator*(), bamHeader);
  //fSVScorer.setSampleCount(scorer2, 1, 1);
  scorer.scoreSV(
      svData,
      mjAssemblyData,
      mjSV,
      mjSVId,
      isJunctionFiltered,
      false,
      true,
      mjModelScoreInfo,
      mjJointModelScoreInfo,
      isMJEvent,
      svEvidenceWriterData);
  BOOST_REQUIRE(isMJEvent);
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
