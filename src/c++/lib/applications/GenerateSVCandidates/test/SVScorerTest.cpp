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
#include "test/testFileMakers.hh"

#include "SVScorer.hh"
#include "SVScorer.cpp"

/// TestSVScorer is a friend of SVScorer. So that can access private
/// method of SVScorer
struct TestSVScorer
{
    void getBreakendMaxMappedDepthAndMQ0(SVScorer& scorer, const bool isTumorOnly, const bool isMaxDepth,
                                         const double cutoffDepth, const SVBreakend& breakend, unsigned& maxDepth,
                                         float& MQ0Frac)
    {
        scorer.getBreakendMaxMappedDepthAndMQ0(isTumorOnly, isMaxDepth, cutoffDepth, breakend, maxDepth, MQ0Frac);
    }

    void computeAllScoreModels(SVScorer& scorer, const bool isSomatic,  const bool isTumorOnly,
                               const std::vector<JunctionCallInfo>& junctionData,
                               SVModelScoreInfo& modelScoreInfo)
    {
        if (isSomatic)
        {
            scorer._sampleCount = 2;
            scorer._diploidSampleCount = 1;
        }
        scorer.computeAllScoreModels(isSomatic, isTumorOnly, junctionData, modelScoreInfo);
    }

    void setSampleCount(SVScorer& scorer, unsigned sampleCount, unsigned diploidSampleCount)
    {
        scorer._sampleCount = sampleCount;
        scorer._diploidSampleCount = diploidSampleCount;
    }
};

BOOST_AUTO_TEST_SUITE( SVScorer_test_suite )

// Create Temporary bam streams of a bam file which contains
// twenty one bam records. Ou of twenty one bam records, 18 are
// with mapping quality 0.
struct BamStream
{
    BamStream()
    {
        const bam_header_info bamHeader(buildTestBamHeader());

        std::string querySeq = "TCTATCACCCATTTTACCACTCACGGGAGCTCTCC";
        bam_record bamRecord1;
        buildTestBamRecord(bamRecord1, 0, 9, 0, 100, 200, 15, "35M", querySeq);
        bamRecord1.set_qname("bamRecord1");
        // small Fragment length read
        bam_record bamRecord2;
        buildTestBamRecord(bamRecord2, 0, 109, 0, 125, 49, 15, "35M", querySeq);
        bamRecord2.set_qname("bamRecord2");
        bam_record bamRecord3;
        buildTestBamRecord(bamRecord3, 0, 109, 0, 180, 100, 15, "35M", querySeq);
        bamRecord3.set_qname("bamRecord3");
        readsToAdd.push_back(bamRecord1);
        readsToAdd.push_back(bamRecord2);
        readsToAdd.push_back(bamRecord3);
        for (unsigned i(0); i < 18; i++)
        {
            bam_record bamRecord;
            buildTestBamRecord(bamRecord, 0, 109, 0, 180, 100, 0, "35M", querySeq);
            std::cout << ((int)bamRecord.map_qual()) << std::endl;
            std::string name = "bamRecord" + std::to_string(i+4);
            bamRecord.set_qname(name.c_str());
            readsToAdd.push_back(bamRecord);
        }
        bamFilename = _bamFilename();
        buildTestBamFile(bamHeader, readsToAdd, bamFilename);

        const std::string referenceFilename = getTestReferenceFilename();
        std::vector<std::string> bamFilenames = { bamFilename };
        std::vector<std::shared_ptr<bam_streamer>> bamStreams;
        openBamStreams(referenceFilename, bamFilenames, bamStreams);
        bamStream = bamStreams[0];
    }

    std::shared_ptr<bam_streamer> bamStream;
    std::vector<bam_record> readsToAdd;
    std::string bamFilename;

private:

    const std::string&
    _bamFilename() const
    {
        return _bamFilenameMaker.getFilename();
    }
    const BamFilenameMaker _bamFilenameMaker;
};

BOOST_FIXTURE_TEST_SUITE( SVScorer_test_suite, BamStream )

// test depth on each location i.e. number of
// read bases overlap in a location.
// We have taken here small fragment legth read as
// here we are testing the utility.
BOOST_AUTO_TEST_CASE( test_addReadToDepthEst )
{
    // read length = 15
    bam_record bamRecord1;
    buildTestBamRecord(bamRecord1, 0, 200, 0, 210, 20, 15, "15M");
    bamRecord1.set_qname("Read-1");

    // Read length = 15
    bam_record bamRecord2;
    buildTestBamRecord(bamRecord2, 0, 210, 0, 220, 20, 15, "15M");
    bamRecord2.set_qname("Read-2");

    std::vector<unsigned>depth(30);
    addReadToDepthEst(bamRecord1, 200, depth);
    addReadToDepthEst(bamRecord2, 200, depth);

    // test the coverage
    for (unsigned i = 0; i < 25; i++)
    {
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
BOOST_AUTO_TEST_CASE( test_lnToProb )
{
    float lower = 5;
    float higher = 5;
    // exp(a-b) = 1. so a = 1/2 and b = 1/2
    lnToProb(lower, higher);
    BOOST_REQUIRE_EQUAL(lower, 0.5);
    BOOST_REQUIRE_EQUAL(higher, 0.5);
}

// Test the following cases
// 1. If forced support is not allowed, either breakpoint-1 or breakpoint-2 should 
//    support split read evidence, otherwise case-2 and case-3 are calculated.
// 2. Alt lnlhood should be maximum of bp1 lnlhood and bp2 lnlhood of alt allele if
//    both the breakpoints support split evidence, otherwise corresponding breakpoint 
//    lnlhood will be taken. 
// 3. Ref lnlhood should be maximum of bp1 lnlhood and bp2 lnlhood of Ref allele if
//    both the breakpoints support split evidence, otherwise corresponding breakpoint 
//    lnlhood will be taken.
BOOST_AUTO_TEST_CASE( test_getSampleSplitReadLnLhood )
{
    float altSplitLnLhoodBP1 = -0.5;
    float altSplitLnLhoodBP2 = -0.8;
    float refSplitLnLhoodBP1 = -0.2;
    float refSplitLnLhoodBP2 = -0.4;

    SVFragmentEvidence fragmentEvidence;

    // BP1 and BP2 of alt are supporting split read evidence where as ref allele
    // is not supporting split read evidence.
    fragmentEvidence.alt.bp1.read1.isSplitSupport = true;
    fragmentEvidence.alt.bp1.read1.splitLnLhood = altSplitLnLhoodBP1;
    fragmentEvidence.alt.bp2.read1.isSplitSupport = true;
    fragmentEvidence.alt.bp2.read1.splitLnLhood = altSplitLnLhoodBP2;
    float refLnLhood;
    float altLnLhood;
    BOOST_REQUIRE(getSampleSplitReadLnLhood(fragmentEvidence, true, refLnLhood, altLnLhood));
    BOOST_REQUIRE_EQUAL(refLnLhood, 0);
    // maximum of (altSplitLnLhoodBP1, altSplitLnLhoodBP2)
    BOOST_REQUIRE_EQUAL(altLnLhood, altSplitLnLhoodBP1);
    
    // BP1 and BP2 of both alt and ref are supporting split read evidence
    fragmentEvidence.ref.bp1.read1.isSplitSupport = true;
    fragmentEvidence.ref.bp1.read1.splitLnLhood = refSplitLnLhoodBP1;
    fragmentEvidence.ref.bp2.read1.isSplitSupport = true;
    fragmentEvidence.ref.bp2.read1.splitLnLhood = refSplitLnLhoodBP2;
    BOOST_REQUIRE(getSampleSplitReadLnLhood(fragmentEvidence, true, refLnLhood, altLnLhood));
    BOOST_REQUIRE_EQUAL(refLnLhood, refSplitLnLhoodBP1);
    BOOST_REQUIRE_EQUAL(altLnLhood, altSplitLnLhoodBP1);
    
    // BP1 and BP2 of ref are supporting split read evidence where as alt allele
    // is not supporting split read evidence.
    fragmentEvidence.alt.bp1.read1.isSplitSupport = false;
    fragmentEvidence.alt.bp1.read1.splitLnLhood = 0;
    fragmentEvidence.alt.bp2.read1.isSplitSupport = false;
    fragmentEvidence.alt.bp2.read1.splitLnLhood = 0;
    BOOST_REQUIRE(getSampleSplitReadLnLhood(fragmentEvidence, true, refLnLhood, altLnLhood));
    BOOST_REQUIRE_EQUAL(refLnLhood, refSplitLnLhoodBP1);
    BOOST_REQUIRE_EQUAL(altLnLhood, 0);
    
    // Both the breakpoints are not supporting split read evidence 
    // for ref and alt alllele. API returns the default value of refLnLhood
    // and altLnLhood which is 1.
    fragmentEvidence.alt.bp1.read1.isSplitSupport = false;
    fragmentEvidence.alt.bp1.read1.splitLnLhood = 0;
    fragmentEvidence.alt.bp2.read1.isSplitSupport = false;
    fragmentEvidence.alt.bp2.read1.splitLnLhood = 0;
    fragmentEvidence.ref.bp1.read1.isSplitSupport = false;
    fragmentEvidence.ref.bp1.read1.splitLnLhood = 0;
    fragmentEvidence.ref.bp2.read1.isSplitSupport = false;
    fragmentEvidence.ref.bp2.read1.splitLnLhood = 0;
    BOOST_REQUIRE(!getSampleSplitReadLnLhood(fragmentEvidence, true, refLnLhood, altLnLhood));
    BOOST_REQUIRE_EQUAL(refLnLhood, 1.);
    BOOST_REQUIRE_EQUAL(altLnLhood, 1.);
}

// Test the following cases:
// 1. If altLnLhood greater than refLnLhood and normalized probability
//    (explained in test_lnToProb) of altLnLhood greater than splitSupportProb(0.999f), 
//    then it is a confident split read count for alt allele.
// 2. if If altLnLhood greater than refLnLhood but normalized probability
//    (explained in test_lnToProb) of altLnLhood is not greater than splitSupportProb(0.999f),
//    then it is not a confident split read count for alt allele. 
// 3. If refLnLhood greater than altLnLhood and normalized probability
//    of altLnLhood greater than splitSupportProb(0.999f), then it is
//    a confident split read count for ref allele.
// 4. For ref allele, if above case-3 is satisfied, then track the confident
//    split read count of BP1 and BP2.
BOOST_AUTO_TEST_CASE( test_addConservativeSplitReadSupport )
{
    // Designed the case-1
    SVSampleInfo sampleInfo1;
    SVFragmentEvidence fragmentEvidence;
    fragmentEvidence.alt.bp1.read1.isSplitSupport = true;
    fragmentEvidence.alt.bp1.read1.splitLnLhood = 6.9077;
    addConservativeSplitReadSupport(fragmentEvidence, true, sampleInfo1);
    BOOST_REQUIRE_EQUAL(sampleInfo1.alt.confidentSplitReadCount, 1);
    
    // Designed the case-2
    SVSampleInfo sampleInfo2;
    fragmentEvidence.alt.bp1.read1.splitLnLhood = 0.2;
    addConservativeSplitReadSupport(fragmentEvidence, true, sampleInfo2);
    BOOST_REQUIRE_EQUAL(sampleInfo2.alt.confidentSplitReadCount, 0);
    
    // Designed the case-3
    // Designed the case-4 for breakpoint-1
    SVSampleInfo sampleInfo3;
    fragmentEvidence.alt.bp1.read1.splitLnLhood = 0;
    fragmentEvidence.ref.bp1.read1.splitLnLhood = 6.9077;
    fragmentEvidence.ref.bp1.read1.isSplitSupport = true;
    addConservativeSplitReadSupport(fragmentEvidence, true, sampleInfo3);
    BOOST_REQUIRE_EQUAL(sampleInfo3.ref.confidentSplitReadCount, 1);
    BOOST_REQUIRE_EQUAL(sampleInfo3.ref.confidentSplitReadAndPairCountRefBp1, 1);
    
    // Designed the case-3
    // Designed the case-4 for breakpoint-2
    SVSampleInfo sampleInfo4;
    fragmentEvidence.alt.bp1.read1.splitLnLhood = 0;
    fragmentEvidence.ref.bp1.read1.splitLnLhood = 0;
    fragmentEvidence.ref.bp1.read1.isSplitSupport = false;
    fragmentEvidence.ref.bp2.read1.splitLnLhood = 6.9077;
    fragmentEvidence.ref.bp2.read1.isSplitSupport = true;
    addConservativeSplitReadSupport(fragmentEvidence, true, sampleInfo4);
    BOOST_REQUIRE_EQUAL(sampleInfo4.ref.confidentSplitReadCount, 1);
    BOOST_REQUIRE_EQUAL(sampleInfo4.ref.confidentSplitReadAndPairCountRefBp2, 1);
}

// Test the following cases:
// 1. If a fragment supports only breakpoint-1, return fragment length
//    probability of breakpoint-1
// 2. If a fragment supports only breakpoint-2, return fragment length
//    probability of breakpoint-2
// 3. If a fragment supports both breakpoint-1 and breakpoin-2, return maximum 
//    of fragment length probability of breakpoint-1 and breakpoint-2
// 4. If a fragment does not support breakpoint-1 or breakpoint-2, return 0. 
BOOST_AUTO_TEST_CASE( test_getSpanningPairAlleleLhood )
{
    float fragLengthProb1 = 0.4;
    float fragLengthProb2 = 0.6;
    
    // Designed the case-1
    SVFragmentEvidenceAllele fragmentEvidenceAllele;
    fragmentEvidenceAllele.bp1.isFragmentSupport = true;
    fragmentEvidenceAllele.bp1.fragLengthProb = fragLengthProb1;
    BOOST_REQUIRE_EQUAL(getSpanningPairAlleleLhood(fragmentEvidenceAllele), fragLengthProb1);
    
    // Designed the case-2
    fragmentEvidenceAllele.bp1.isFragmentSupport = false;
    fragmentEvidenceAllele.bp2.isFragmentSupport = true;
    fragmentEvidenceAllele.bp2.fragLengthProb = fragLengthProb2;
    BOOST_REQUIRE_EQUAL(getSpanningPairAlleleLhood(fragmentEvidenceAllele), fragLengthProb2);

    // Designed the case-3
    fragmentEvidenceAllele.bp1.isFragmentSupport = true;
    BOOST_REQUIRE_EQUAL(getSpanningPairAlleleLhood(fragmentEvidenceAllele), fragLengthProb2);

    // Designed the case-4
    fragmentEvidenceAllele.bp1.isFragmentSupport = false;
    fragmentEvidenceAllele.bp2.isFragmentSupport = false;
    BOOST_REQUIRE_EQUAL(getSpanningPairAlleleLhood(fragmentEvidenceAllele), 0);
}

// Test Spanning pair count based on whether a fragment supports
// breakpoint-1(BP1) or breakpoint-2(BP2) for alt and ref allele.
BOOST_AUTO_TEST_CASE( test_addSpanningPairSupport )
{
    // Fragment does not support BP1 or BP2 for both alt and ref allele
    SVFragmentEvidence fragmentEvidence;
    SVSampleInfo sampleInfo1;
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
// 2. When fraction of length probability of alt allele greater than fraction of 
//    length probability of ref allele. confidentSemiMappedSpanningPairCount of alt allele is
//    incremented by 1.
// 3. When fully mapped fragment means both read1 and read2 both are anchored read. So 
//    confidentSpanningPairCount is incremented by 1.
// 4. When fraction of length probability of ref allele greater than fraction of length 
//    probability of alt allele. So confidentSemiMappedSpanningPairCount of ref allele is
//    incremented by 1.
// 5. when fraction of length probability of ref allele greater than alt allele and
//    fragment suypports BP1 or BP2 of ref allele. Either confidentSplitReadAndPairCountRefBp1 or 
//    confidentSplitReadAndPairCountRefBp2 or both is(are) incremented by 1. 
// 6. When fragment length probability = 0. It throws an exception.
BOOST_AUTO_TEST_CASE( test_addConservativeSpanningPairSupport )
{
    float fragLengthProb1 = 0.4;
    // Designed the case-1
    SVFragmentEvidence fragmentEvidence1;
    SVSampleInfo sampleInfo1;
    addConservativeSpanningPairSupport(fragmentEvidence1, sampleInfo1);
    BOOST_REQUIRE_EQUAL(sampleInfo1.alt.confidentSemiMappedSpanningPairCount, 0);
    BOOST_REQUIRE_EQUAL(sampleInfo1.alt.confidentSpanningPairCount, 0);
    BOOST_REQUIRE_EQUAL(sampleInfo1.ref.confidentSemiMappedSpanningPairCount, 0);
    BOOST_REQUIRE_EQUAL(sampleInfo1.ref.confidentSpanningPairCount, 0);
    BOOST_REQUIRE_EQUAL(sampleInfo1.ref.confidentSplitReadAndPairCountRefBp1, 0);
    BOOST_REQUIRE_EQUAL(sampleInfo1.ref.confidentSplitReadAndPairCountRefBp2, 0);
    
    // Designed the case-2
    SVSampleInfo sampleInfo2;
    fragmentEvidence1.alt.bp1.isFragmentSupport = true;
    fragmentEvidence1.alt.bp1.fragLengthProb = fragLengthProb1;
    addConservativeSpanningPairSupport(fragmentEvidence1, sampleInfo2);
    BOOST_REQUIRE_EQUAL(sampleInfo2.alt.confidentSemiMappedSpanningPairCount, 1);
    BOOST_REQUIRE_EQUAL(sampleInfo2.alt.confidentSpanningPairCount, 0);
    BOOST_REQUIRE_EQUAL(sampleInfo2.ref.confidentSemiMappedSpanningPairCount, 0);
    BOOST_REQUIRE_EQUAL(sampleInfo2.ref.confidentSpanningPairCount, 0);
    BOOST_REQUIRE_EQUAL(sampleInfo2.ref.confidentSplitReadAndPairCountRefBp1, 0);
    BOOST_REQUIRE_EQUAL(sampleInfo2.ref.confidentSplitReadAndPairCountRefBp2, 0);
    
    // Designed the case-3
    SVSampleInfo sampleInfo3;
    fragmentEvidence1.alt.bp1.isFragmentSupport = true;
    fragmentEvidence1.alt.bp1.fragLengthProb = fragLengthProb1;
    fragmentEvidence1.read1.isScanned = true;
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
    SVFragmentEvidence fragmentEvidence2;
    SVSampleInfo sampleInfo4;
    fragmentEvidence2.ref.bp1.isFragmentSupport = true;
    fragmentEvidence2.ref.bp1.fragLengthProb = fragLengthProb1;
    addConservativeSpanningPairSupport(fragmentEvidence2, sampleInfo4);
    BOOST_REQUIRE_EQUAL(sampleInfo4.alt.confidentSemiMappedSpanningPairCount, 0);
    BOOST_REQUIRE_EQUAL(sampleInfo4.alt.confidentSpanningPairCount, 0);
    BOOST_REQUIRE_EQUAL(sampleInfo4.ref.confidentSemiMappedSpanningPairCount, 1);
    BOOST_REQUIRE_EQUAL(sampleInfo4.ref.confidentSpanningPairCount, 0);
    BOOST_REQUIRE_EQUAL(sampleInfo4.ref.confidentSplitReadAndPairCountRefBp1, 0);
    BOOST_REQUIRE_EQUAL(sampleInfo4.ref.confidentSplitReadAndPairCountRefBp2, 0);
    
    // Designed the case-3 (fully mapped fragment) and the case-5 when 
    // fragment supports BP1 of ref allele
    SVSampleInfo sampleInfo5;
    fragmentEvidence2.ref.bp1.isFragmentSupport = true;
    fragmentEvidence2.ref.bp1.fragLengthProb = fragLengthProb1;
    fragmentEvidence2.read1.isScanned = true;
    fragmentEvidence2.read1.setAnchored(true);
    fragmentEvidence2.read2.isScanned = true;
    fragmentEvidence2.read2.setAnchored(true);
    addConservativeSpanningPairSupport(fragmentEvidence2, sampleInfo5);
    BOOST_REQUIRE_EQUAL(sampleInfo5.alt.confidentSemiMappedSpanningPairCount, 0);
    BOOST_REQUIRE_EQUAL(sampleInfo5.alt.confidentSpanningPairCount, 0);
    BOOST_REQUIRE_EQUAL(sampleInfo5.ref.confidentSemiMappedSpanningPairCount, 1);
    BOOST_REQUIRE_EQUAL(sampleInfo5.ref.confidentSpanningPairCount, 1);
    BOOST_REQUIRE_EQUAL(sampleInfo5.ref.confidentSplitReadAndPairCountRefBp1, 1);
    BOOST_REQUIRE_EQUAL(sampleInfo5.ref.confidentSplitReadAndPairCountRefBp2, 0);
    
    // Designed the case-3 (fully mapped fragment) and the case-5 when 
    // fragment supports BP2 of ref allele
    SVSampleInfo sampleInfo6;
    fragmentEvidence2.ref.bp1.isFragmentSupport = false;
    fragmentEvidence2.ref.bp2.isFragmentSupport = true;
    fragmentEvidence2.ref.bp2.fragLengthProb = fragLengthProb1;
    fragmentEvidence2.read1.isScanned = true;
    fragmentEvidence2.read1.setAnchored(true);
    fragmentEvidence2.read2.isScanned = true;
    fragmentEvidence2.read2.setAnchored(true);
    addConservativeSpanningPairSupport(fragmentEvidence2, sampleInfo6);
    BOOST_REQUIRE_EQUAL(sampleInfo6.alt.confidentSemiMappedSpanningPairCount, 0);
    BOOST_REQUIRE_EQUAL(sampleInfo6.alt.confidentSpanningPairCount, 0);
    BOOST_REQUIRE_EQUAL(sampleInfo6.ref.confidentSemiMappedSpanningPairCount, 1);
    BOOST_REQUIRE_EQUAL(sampleInfo6.ref.confidentSpanningPairCount, 1);
    BOOST_REQUIRE_EQUAL(sampleInfo6.ref.confidentSplitReadAndPairCountRefBp1, 0);
    BOOST_REQUIRE_EQUAL(sampleInfo6.ref.confidentSplitReadAndPairCountRefBp2, 1);
    
    // Designed the case-6. It will throw an exception
    fragmentEvidence2.ref.bp2.fragLengthProb = 0;
    BOOST_CHECK_THROW(addConservativeSpanningPairSupport(fragmentEvidence2, sampleInfo6), illumina::common::GeneralException);
}

// Test test_addConservativeSpanningPairSupport and test_addConservativeSplitReadSupport
// for a sample means api will track all the stats of a sample.
BOOST_AUTO_TEST_CASE( test_getSampleCounts )
{
    // sample evidence information
    SVEvidence::evidenceTrack_t samples;
    SVFragmentEvidence fragmentEvidence1;
    SVSampleInfo sampleInfo;
    fragmentEvidence1.alt.bp1.isFragmentSupport = true;
    fragmentEvidence1.alt.bp1.fragLengthProb = 0.4f;
    fragmentEvidence1.read1.isScanned = true;
    fragmentEvidence1.read1.setAnchored(true);
    fragmentEvidence1.read2.isScanned = true;
    fragmentEvidence1.read2.setAnchored(true);
    SVEvidence evidence;
    samples["Fragment-1"] = fragmentEvidence1;
    SVFragmentEvidence fragmentEvidence2;
    fragmentEvidence2.alt.bp1.read1.isSplitSupport = true;
    fragmentEvidence2.alt.bp1.read1.splitLnLhood = 6.9077;
    samples["Fragment-2"] = fragmentEvidence2;
    getSampleCounts(samples, sampleInfo);
    
    // following two cases are coming from fragment-1 in test_addConservativeSpanningPairSupport
    BOOST_REQUIRE_EQUAL(sampleInfo.alt.confidentSpanningPairCount, 1);
    BOOST_REQUIRE_EQUAL(sampleInfo.alt.confidentSemiMappedSpanningPairCount, 1);
    // This one is coming from fragment-2 in test_addConservativeSplitReadSupport
    BOOST_REQUIRE_EQUAL(sampleInfo.alt.confidentSplitReadCount, 1);
}

// Test test_getSampleCounts for multiple samples
BOOST_AUTO_TEST_CASE( test_getSVSupportSummary )
{
    // Data for sample-1
    SVEvidence::evidenceTrack_t sample1;
    SVFragmentEvidence fragmentEvidence1;
    fragmentEvidence1.alt.bp1.isFragmentSupport = true;
    fragmentEvidence1.alt.bp1.fragLengthProb = 0.4f;
    fragmentEvidence1.read1.isScanned = true;
    fragmentEvidence1.read1.setAnchored(true);
    fragmentEvidence1.read2.isScanned = true;
    fragmentEvidence1.read2.setAnchored(true);
    sample1["Fragment-1"] = fragmentEvidence1;
    SVFragmentEvidence fragmentEvidence2;
    fragmentEvidence2.alt.bp1.read1.isSplitSupport = true;
    fragmentEvidence2.alt.bp1.read1.splitLnLhood = 6.9077;
    sample1["Fragment-2"] = fragmentEvidence2;
    SVEvidence evidence;
    evidence.samples.push_back(sample1);
    
    // Data for sample-2
    SVEvidence::evidenceTrack_t sample2;
    SVFragmentEvidence fragmentEvidence3;
    fragmentEvidence3.ref.bp1.isFragmentSupport = true;
    fragmentEvidence3.ref.bp1.fragLengthProb = 0.4f;
    fragmentEvidence3.read1.isScanned = true;
    fragmentEvidence3.read1.setAnchored(true);
    fragmentEvidence3.read2.isScanned = true;
    fragmentEvidence3.read2.setAnchored(true);
    sample2["Fragment-3"] = fragmentEvidence3;
    SVFragmentEvidence fragmentEvidence4;
    fragmentEvidence4.ref.bp1.read1.isSplitSupport = true;
    fragmentEvidence4.ref.bp1.read1.splitLnLhood = 6.9077;
    sample2["Fragment-4"] = fragmentEvidence4;
    evidence.samples.push_back(sample2);

    SVScoreInfo scoreInfo;
    scoreInfo.samples.resize(2);
    getSVSupportSummary(evidence, scoreInfo);
    
    // check result for sample-1 with respect to test_getSampleCounts
    BOOST_REQUIRE_EQUAL(scoreInfo.samples[0].alt.confidentSpanningPairCount, 1);
    BOOST_REQUIRE_EQUAL(scoreInfo.samples[0].alt.confidentSemiMappedSpanningPairCount, 1);
    BOOST_REQUIRE_EQUAL(scoreInfo.samples[0].alt.confidentSplitReadCount, 1);
    // check result for sample-2 with respect to test_getSampleCounts
    BOOST_REQUIRE_EQUAL(scoreInfo.samples[1].ref.confidentSpanningPairCount, 1);
    BOOST_REQUIRE_EQUAL(scoreInfo.samples[1].ref.confidentSemiMappedSpanningPairCount, 1);
    BOOST_REQUIRE_EQUAL(scoreInfo.samples[1].ref.confidentSplitReadCount, 1);
}

// if there's a difference in fragment support for one "frag-favored" allele, then
// there must also be either neutral split support or split support in favor of the same allele.
// So we can reset the all values for this fragment. Following cases need to be tested:
// 1. When refSplitLnLhood > altSplitLnLhood at read1, but altPairLhood > refPairLhood
// 2. When altSplitLnLhood > refSplitLnLhood at read1, but refPairLhood > altPairLhood
// 3. When refSplitLnLhood > altSplitLnLhood at read2, but altPairLhood > refPairLhood
// 4. When altSplitLnLhood > refSplitLnLhood at read2, but refPairLhood > altPairLhood
BOOST_AUTO_TEST_CASE( test_resolvePairSplitConflicts )
{
    float fragLengthProb1 = 0.4;
    SVCandidate candidate;
    candidate.setPrecise();
    SVEvidence::evidenceTrack_t evidenceTrack;
    std::string fragLabel = "frag-1";
    SVFragmentEvidence fragmentEvidence1;
    fragmentEvidence1.alt.bp1.isFragmentSupport = true;
    fragmentEvidence1.alt.bp1.fragLengthProb = fragLengthProb1;
    fragmentEvidence1.ref.bp1.read1.isSplitSupport = true;
    fragmentEvidence1.ref.bp1.read1.splitLnLhood = 0.2;
    evidenceTrack[fragLabel] = fragmentEvidence1;
    SVEvidence evidence;
    evidence.samples.resize(1);
    evidence.samples.push_back(evidenceTrack);
    // Designed the case-1
    resolvePairSplitConflicts(candidate, evidence);
    BOOST_REQUIRE(!evidence.getSampleEvidence(0)[fragLabel].alt.bp1.isFragmentSupport);
    BOOST_REQUIRE_EQUAL(evidence.getSampleEvidence(0)[fragLabel].alt.bp1.fragLengthProb, 0);

    SVFragmentEvidence fragmentEvidence2;
    fragmentEvidence2.ref.bp1.isFragmentSupport = true;
    fragmentEvidence2.ref.bp1.fragLengthProb = fragLengthProb1;
    fragmentEvidence2.alt.bp1.read1.isSplitSupport = true;
    fragmentEvidence2.alt.bp1.read1.splitLnLhood = 0.2;
    evidenceTrack[fragLabel] = fragmentEvidence2;
    evidence.samples.clear();
    evidence.samples.push_back(evidenceTrack);
    // Designed the case-2
    resolvePairSplitConflicts(candidate, evidence);
    BOOST_REQUIRE(!evidence.getSampleEvidence(0)[fragLabel].ref.bp1.isFragmentSupport);
    BOOST_REQUIRE_EQUAL(evidence.getSampleEvidence(0)[fragLabel].ref.bp1.fragLengthProb, 0);
    
    // Designed the case-3
    SVFragmentEvidence fragmentEvidence3;
    fragmentEvidence3.ref.bp1.isFragmentSupport = true;
    fragmentEvidence3.ref.bp1.fragLengthProb = fragLengthProb1;
    fragmentEvidence3.alt.bp1.read2.isSplitSupport = true;
    fragmentEvidence3.alt.bp1.read2.splitLnLhood = 0.2;
    evidenceTrack[fragLabel] = fragmentEvidence3;
    evidence.samples.clear();
    evidence.samples.push_back(evidenceTrack);
    resolvePairSplitConflicts(candidate, evidence);
    BOOST_REQUIRE(!evidence.getSampleEvidence(0)[fragLabel].ref.bp1.isFragmentSupport);
    BOOST_REQUIRE_EQUAL(evidence.getSampleEvidence(0)[fragLabel].ref.bp1.fragLengthProb, 0);
    
    // Designed the case-4
    SVFragmentEvidence fragmentEvidence4;
    fragmentEvidence4.alt.bp1.isFragmentSupport = true;
    fragmentEvidence4.alt.bp1.fragLengthProb = fragLengthProb1;
    fragmentEvidence4.ref.bp1.read2.isSplitSupport = true;
    fragmentEvidence4.ref.bp1.read2.splitLnLhood = 0.2;
    evidenceTrack[fragLabel] = fragmentEvidence4;
    evidence.samples.clear();
    evidence.samples.push_back(evidenceTrack);
    resolvePairSplitConflicts(candidate, evidence);
    BOOST_REQUIRE(!evidence.getSampleEvidence(0)[fragLabel].ref.bp1.isFragmentSupport);
    BOOST_REQUIRE_EQUAL(evidence.getSampleEvidence(0)[fragLabel].ref.bp1.fragLengthProb, 0);
}

// Test the following utility formula:
// log(selfChimeraProb.comp*fragProb + otherChimeraProb.prob)*power,
// where fragProb is mention in test_getSpanningPairAlleleLhood and 
// power = 2.
BOOST_AUTO_TEST_CASE( test_incrementSpanningPairAlleleLnLhood )
{
    float fragLengthProb1 = 0.4;
    float fragLengthProb2 = 0.6;
    ProbSet selfChimeraProb(0.5);
    ProbSet otherChimeraProb(0.2);
    // check for breakpoint-1
    SVFragmentEvidenceAllele fragmentEvidenceAllele;
    fragmentEvidenceAllele.bp1.isFragmentSupport = true;
    fragmentEvidenceAllele.bp1.fragLengthProb = fragLengthProb1;
    double bpLnLhood1Expected = -1.83258;
    double  bpLnLhood1 = 0;
    incrementSpanningPairAlleleLnLhood(selfChimeraProb, otherChimeraProb, fragmentEvidenceAllele, 2, bpLnLhood1);
    // Checking upto 5 decimal places
    BOOST_REQUIRE_EQUAL(((int)(bpLnLhood1 * 100000))/ 100000.0, bpLnLhood1Expected);
    
    // check for breakpoint-2
    double bpLnLhood2Expected = -1.38629;
    double  bpLnLhood2 = 0;
    fragmentEvidenceAllele.bp1.isFragmentSupport = false;
    fragmentEvidenceAllele.bp2.isFragmentSupport = true;
    fragmentEvidenceAllele.bp2.fragLengthProb = fragLengthProb2;
    incrementSpanningPairAlleleLnLhood(selfChimeraProb, otherChimeraProb, fragmentEvidenceAllele, 2, bpLnLhood2);
    // checking upto 5 decimal places
    BOOST_REQUIRE_EQUAL(((int)(bpLnLhood2 * 100000))/ 100000.0, bpLnLhood2Expected);
}

// fragLnLHood = log_sum(selfMapProb.lnComp+alignLnLhood), otherMapProb.lnProb),
//     log_sum(x1,x2) = x1 + log(1+exp(x2-x1)) if x2>x1
//                    = x2 + log(1+exp(x1-x2)) if x1>x2
//     lnComp = log(1 - lnProb)
// Following cases need to tested
// 1. When a fragment does not support none of the breakpoints
// 2. When a fragment support both the breakpoints
// 3. When a fragment support both the breakpoints and splitLnLhood of BP1
//    greater than splitLnLhood of BP2
// 4. When a fragment support both the breakpoints and splitLnLhood of BP2
//    greater than splitLnLhood of BP1
BOOST_AUTO_TEST_CASE( test_incrementAlleleSplitReadLhood )
{
    ProbSet selfChimeraProb(0.5);
    ProbSet otherChimeraProb(0.2);
    // Designed the case-1
    SVFragmentEvidenceAllele fragmentEvidenceAllele;
    fragmentEvidenceAllele.bp1.read1.splitLnLhood = 6;
    fragmentEvidenceAllele.bp2.read1.splitLnLhood = 5;
    bool isReadEvaluated;
    double expectedValue1 = ((int)(4.30954 * 100000)) / 100000.0;
    double expectedValue2 = ((int)(5.30784 * 100000)) / 100000.0;
    BOOST_REQUIRE_EQUAL(((int)(incrementAlleleSplitReadLhood(selfChimeraProb, otherChimeraProb,
                                                             fragmentEvidenceAllele, 0.5,
                                                             fragmentEvidenceAllele.isAnySplitReadSupport(true),
                                                             true, isReadEvaluated)*100000))/100000.0, expectedValue1);

    // Designed the case-2 and the case-3 when splitLnLhood of BP1
    // greater than splitLnLhood of BP2. Here alignLnLhood = 2
    fragmentEvidenceAllele.bp1.read1.isSplitSupport = true;
    fragmentEvidenceAllele.bp2.read1.isSplitSupport = true;
    BOOST_REQUIRE_EQUAL(((int)(incrementAlleleSplitReadLhood(selfChimeraProb, otherChimeraProb,
                                                             fragmentEvidenceAllele, 0.5,
                                                             fragmentEvidenceAllele.isAnySplitReadSupport(true),
                                                             true, isReadEvaluated)*100000))/100000.0, expectedValue2);
    // Designed case-2 and case-3 when splitLnLhood of BP2
    // greater than splitLnLhood of BP1. Here alignLnLhood = 5
    fragmentEvidenceAllele.bp1.read1.splitLnLhood = 2;
    fragmentEvidenceAllele.bp2.read1.splitLnLhood = 5;
    BOOST_REQUIRE_EQUAL(((int)(incrementAlleleSplitReadLhood(selfChimeraProb, otherChimeraProb,
                                                             fragmentEvidenceAllele, 0.5,
                                                             fragmentEvidenceAllele.isAnySplitReadSupport(true),
                                                             true, isReadEvaluated)*100000))/100000.0, expectedValue1);
}

// Following cases need to be tested
// 1. When a fragment does not support none of the breakpoints
// 2. When a fragment does not support any tier2 of none of the breakpoints
// 3. When a fragment supports one of the breakpoints
BOOST_AUTO_TEST_CASE( test_incrementSplitReadLhood )
{
    ProbSet refMapProb(0.5);
    ProbSet altMapProb(0.2);
    
    // Designed the case-1. In this case refSplitLnLhood and altSplitLnLhood are 0 as
    // fragment does not support any of the breakpoints
    SVFragmentEvidence fragmentEvidence;
    fragmentEvidence.alt.bp1.read1.splitLnLhood = 6;
    fragmentEvidence.alt.bp2.read1.splitLnLhood = 5;
    bool isReadEvaluated;
    std::string fragLabel = "Frag-1";
    double refSplitLnLhood = 0;
    double altSplitLnLhood = 0;
    incrementSplitReadLhood(fragLabel, fragmentEvidence, refMapProb, altMapProb, false, true, refSplitLnLhood, altSplitLnLhood, isReadEvaluated);
    BOOST_REQUIRE_EQUAL(refSplitLnLhood, 0);
    BOOST_REQUIRE_EQUAL(altSplitLnLhood, 0);
    
    // Designed the case-2 when a fragment does not support any tier2
    incrementSplitReadLhood(fragLabel, fragmentEvidence, refMapProb, altMapProb, true, true, refSplitLnLhood, altSplitLnLhood, isReadEvaluated);
    BOOST_REQUIRE_EQUAL(refSplitLnLhood, 0);
    BOOST_REQUIRE_EQUAL(altSplitLnLhood, 0);
    
    // Designed the case-3 which is explained in test_incrementAlleleSplitReadLhood
    fragmentEvidence.alt.bp1.read1.isSplitSupport = true;
    fragmentEvidence.alt.bp2.read1.isSplitSupport = true;
    incrementSplitReadLhood(fragLabel, fragmentEvidence, refMapProb, altMapProb, false, true, refSplitLnLhood, altSplitLnLhood, isReadEvaluated);
    BOOST_REQUIRE(refSplitLnLhood != 0);
    BOOST_REQUIRE(altSplitLnLhood != 0);
}

// Test the following cases:
// 1. When evaluation is done only for read-1. Then api should return sum of value of read-1 and fragPair.
// 2. When evaluation is done only for read-2. Then api should return sum of value of read-2 and fragPair.
// 3. When evaluation is done based on both read-1 and read-2. Then api should return sum of (maximum value of read-1 and 
//    read-2) and fragPair.
// 4. When evaluation is done without read-1 and read-2. Then api should return value of fragPair.
BOOST_AUTO_TEST_CASE( test_getFragLnLhood )
{
    // book keeping of lnLhood of read-1 and read-2
    AlleleLnLhood lnLhood;
    lnLhood.read1Split = 0.2;
    lnLhood.read2Split = 0.5;
    lnLhood.fragPair = 0.4;
    
    // Designed the case-1
    BOOST_REQUIRE_EQUAL(getFragLnLhood(lnLhood, true, false), lnLhood.read1Split + lnLhood.fragPair);
    
    // Designed the case-2
    BOOST_REQUIRE_EQUAL(getFragLnLhood(lnLhood, false, true), lnLhood.read2Split + lnLhood.fragPair);
    
    // Designed the case-3
    BOOST_REQUIRE_EQUAL(getFragLnLhood(lnLhood, true, true), lnLhood.read2Split + lnLhood.fragPair);
    
    // Designed the case-4
    BOOST_REQUIRE_EQUAL(getFragLnLhood(lnLhood, false, false), lnLhood.fragPair);
}

// Spanning weight is calculated as,
// spanning weight = min(1, max(0, (val-_min)*_factor)),
// where, _min = 300, if sv insertion is not large (<100)
//             = 100, if sv insertion is large(> 100)
//        val = center size of sv, if sv insertion is not large (<100)
//            = insertion size, if sv insertion is large (>100)
//        factor = 1/(_max-_min), 
//        _max = 500, if sv insertion is not large (<100)
//             = 150, if sv insertion is large(> 100)
// Following cases need to tested:
// 1. When sv type is not indel. In this case it will return 1.
// 2. When sv insertion is not large
// 3  when sv insertion is large
BOOST_AUTO_TEST_CASE( test_getSpanningPairWeight )
{
    SVCandidate candidate;
    // Designed the case-1
    candidate.bp1.interval = GenomeInterval(0, 1000, 1200);
    candidate.bp2.interval = GenomeInterval(1, 1000, 1200);
    BOOST_REQUIRE_EQUAL(getSpanningPairWeight(candidate), 1);
    
    // Designed the case-2
    candidate.bp1.interval = GenomeInterval(0, 1000, 1200);
    candidate.bp2.interval = GenomeInterval(0, 500, 600);
    candidate.bp1.state = SVBreakendState::LEFT_OPEN;
    candidate.bp2.state = SVBreakendState::RIGHT_OPEN;
    BOOST_REQUIRE_EQUAL(getSpanningPairWeight(candidate), 1);
    
    // Designed the case-3
    candidate.bp1.interval = GenomeInterval(0, 1000, 1200);
    candidate.bp2.interval = GenomeInterval(0, 980, 990);
    candidate.insertSeq = "GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGT"
                          "ATTTTCGTCTGGGGGGTGTGCACGCGATAGCATTGCGAGACGCTGGA";
    BOOST_REQUIRE_EQUAL(getSpanningPairWeight(candidate), 0.04f);
}

// large SV prior weight is same as explained in test_getSpanningPairWeight. But
// here _min = 5000 and _max = 10000.
// For interchrosomal SV event it will return 1.
BOOST_AUTO_TEST_CASE( test_largeNoiseSVPriorWeight )
{
    SVCandidate candidate;
    // Interchrosomal SV. So largeNoiseSVPriorWeight = 1.
    candidate.bp1.interval = GenomeInterval(0, 1000, 1200);
    candidate.bp2.interval = GenomeInterval(1, 1000, 1200);
    BOOST_REQUIRE_EQUAL(largeNoiseSVPriorWeight(candidate), 1);
    
    // same as explained in test_getSpanningPairWeight
    candidate.bp1.interval = GenomeInterval(0, 1000, 1200);
    candidate.bp2.interval = GenomeInterval(0, 500, 600);
    BOOST_REQUIRE_EQUAL(largeNoiseSVPriorWeight(candidate), 0);
}

// Test the following cases
// 1. return true if any evidence exists for fragment
// 2. LnLhood values are computed as explained in test_incrementSpanningPairAlleleLnLhood
BOOST_AUTO_TEST_CASE( test_getRefAltFromFrag )
{
    float spanningPairWeight = 0.2;
    double semiMappedPower = 2;
    ProbSet refChimeraProb(0.2);
    ProbSet altChimeraProb(0.3);
    ProbSet refSplitMapProb(0.1);
    ProbSet altSplitMapProb(0.5);
    std::string fragLabel = "frag-1";
    SVFragmentEvidence evidence;
    AlleleLnLhood refLnLhoodSet;
    AlleleLnLhood altLnLhoodSet;
    bool isRead1Evaluated;
    bool isRead2Evaluated;
    // evidence is not true for the fragment
    BOOST_REQUIRE(!getRefAltFromFrag(spanningPairWeight, semiMappedPower, refChimeraProb, altChimeraProb, refSplitMapProb, altSplitMapProb,
                                     true, fragLabel, evidence, refLnLhoodSet,altLnLhoodSet, isRead1Evaluated, isRead2Evaluated));
    
    // evidence is true for the fragment
    evidence.read1.isScanned = true;
    evidence.read2.isScanned = true;
    evidence.read1.setAnchored(true);
    evidence.read2.setAnchored(true);
    evidence.alt.bp1.isFragmentSupport = true;
    BOOST_REQUIRE(getRefAltFromFrag(spanningPairWeight, semiMappedPower, refChimeraProb, altChimeraProb, refSplitMapProb, altSplitMapProb,
                                     false, fragLabel, evidence, refLnLhoodSet,altLnLhoodSet, isRead1Evaluated, isRead2Evaluated));
    
    // For semi mapped read.
    evidence.read2.setAnchored(false);
    BOOST_REQUIRE(getRefAltFromFrag(spanningPairWeight, semiMappedPower, refChimeraProb, altChimeraProb, refSplitMapProb, altSplitMapProb,
                                    false, fragLabel, evidence, refLnLhoodSet,altLnLhoodSet, isRead1Evaluated, isRead2Evaluated));
    // check test_incrementSpanningPairAlleleLnLhood
    BOOST_REQUIRE(refLnLhoodSet.fragPair != 0);
    BOOST_REQUIRE(altLnLhoodSet.fragPair != 0);

    // fragLengthProb of alt allle is more than fragLengthProb of ref allele
    evidence.alt.bp1.fragLengthProb = 0.3;
    BOOST_REQUIRE(getRefAltFromFrag(spanningPairWeight, semiMappedPower, refChimeraProb, altChimeraProb, refSplitMapProb, altSplitMapProb,
                                    false, fragLabel, evidence, refLnLhoodSet,altLnLhoodSet, isRead1Evaluated, isRead2Evaluated));
    // check test_incrementSpanningPairAlleleLnLhood
    BOOST_REQUIRE(refLnLhoodSet.fragPair != 0);
    BOOST_REQUIRE(altLnLhoodSet.fragPair != 0);
}

// Test the Diploid scoring model
// For reference and alternate alleles A = {r, x},
// the genotype states at each allele are restricted to G = {rr, rx, xx}.
// So if refLnFragLhood of ref allele is R and altLnFragLhood of alt allele
// is L, then score(S) of a genotype (gt = REF or HET or HOM) is
// S = log_sum(R + lnPreComp(gt), L + LnPre(gt)), where,
//     LnPre(gt) = log(Prior probability of expected alt allele)
//     lnPreComp(gt) = log(1-Prior probability of expected alt allele)
//     log_sum(x1,x2) = x1 + log(1+exp(x2-x1)) if x2>x1
//                    = x2 + log(1+exp(x1-x2)) if x1>x2
BOOST_AUTO_TEST_CASE( test_addDiploidLoglhood )
{
    SVEvidence::evidenceTrack_t evidenceTrack;
    std::string fragLabel = "frag-1";
    SVFragmentEvidence evidence;
    evidenceTrack[fragLabel] = evidence;
    std::array<double,DIPLOID_GT::SIZE> loglhood;
    for (unsigned gt(0); gt<DIPLOID_GT::SIZE; ++gt)
    {
        loglhood[gt] = 0;
    }
    //fragment was not evaluated for pair or split support for either allele
    addDiploidLoglhood(0.2, evidenceTrack, loglhood);
    for (unsigned gt(0); gt<DIPLOID_GT::SIZE; ++gt)
    {
        BOOST_REQUIRE(loglhood[gt] == 0);
    }
    // fragment was evaluated for alt allele
    evidence.read1.isScanned = true;
    evidence.read2.isScanned = true;
    evidence.read1.setAnchored(true);
    evidence.read2.setAnchored(true);
    evidence.alt.bp1.isFragmentSupport = true;
    evidenceTrack[fragLabel] = evidence;
    addDiploidLoglhood(0.2, evidenceTrack, loglhood);
    // above mentioned formula is applied in loglhood
    for (unsigned gt(0); gt<DIPLOID_GT::SIZE; ++gt)
    {
        BOOST_REQUIRE(loglhood[gt] != 0); // just checking not equal to zero because of floating point issue
    }
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
BOOST_AUTO_TEST_CASE( test_DiploidSVLabel )
{
    CallOptionsDiploid diploidOpt;
    const bam_header_info bamHeader(buildTestBamHeader());
    std::unique_ptr<SVLocusScanner> scanner(buildTestSVLocusScanner(bamHeader));
    CallOptionsDiploidDeriv diploidDeriv(diploidOpt);
    // dummy depth file that means depth file is not present
    ChromDepthFilterUtil dFilter("", 2, bamHeader);
    
    // Diploid sample info
    SVScoreInfoDiploid diploidInfo;
    diploidInfo.setSampleCount(1); // diploid sample count = 1
    SVFragmentEvidence fragmentEvidence;
    fragmentEvidence.read1.isScanned = true;
    fragmentEvidence.read2.isScanned = true;
    fragmentEvidence.read1.setAnchored(true);
    fragmentEvidence.read2.setAnchored(true);
    fragmentEvidence.alt.bp1.isFragmentSupport = true;
    SVEvidence::evidenceTrack_t evidenceTrack;
    std::string fragLabel = "frag-1";
    evidenceTrack[fragLabel] = fragmentEvidence;
    SVEvidence evidence;
    evidence.samples.push_back(evidenceTrack);
    
    // Diploid score info
    SVScoreInfo scoreInfo;
    scoreInfo.samples.resize(1);
    SVSampleInfo diploidSampleInfo;
    diploidSampleInfo.alt.confidentSpanningPairCount = 100;
    diploidSampleInfo.alt.confidentSemiMappedSpanningPairCount = 50;
    diploidSampleInfo.alt.confidentSplitReadCount = 10;
    scoreInfo.samples.push_back(diploidSampleInfo);
    SVCandidate candidate;
    JunctionCallInfo callInfo;
    callInfo.init(candidate, evidence, scoreInfo, 0.2);
    std::vector<JunctionCallInfo> junctionData = {callInfo};
    
    scoreDiploidSV(diploidOpt, scanner.operator*(), diploidDeriv, dFilter, junctionData, diploidInfo);
    // Designed the case-1 and case-2
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
    BOOST_REQUIRE_EQUAL(diploidInfo.filters.size(), 3);
    BOOST_REQUIRE_EQUAL(*(diploidInfo.filters.begin()), diploidOpt.maxMQ0FracLabel);
    BOOST_REQUIRE_EQUAL(*(std::next(diploidInfo.filters.begin(), 1)), diploidOpt.minAltFilterLabel);
    BOOST_REQUIRE_EQUAL(*(std::next(diploidInfo.filters.begin(), 2)), diploidOpt.failedSampleFTLabel);

    // Designed the case-1, case-2 and case-4
    candidate.bp1.interval = GenomeInterval(0, 100, 200);
    candidate.bp2.interval = GenomeInterval(1, 200, 300);
    junctionData.clear();
    diploidInfo.filters.clear();
    callInfo.init(candidate, evidence, scoreInfo, 0.2);
    junctionData.push_back(callInfo);
    scoreDiploidSV(diploidOpt, scanner.operator*(), diploidDeriv, dFilter, junctionData, diploidInfo);
    BOOST_REQUIRE_EQUAL(diploidInfo.filters.size(), 3);
    BOOST_REQUIRE_EQUAL(*(diploidInfo.filters.begin()), diploidOpt.minAltFilterLabel);
    BOOST_REQUIRE_EQUAL(*(std::next(diploidInfo.filters.begin(), 1)), diploidOpt.noPairSupportLabel);
    BOOST_REQUIRE_EQUAL(*(std::next(diploidInfo.filters.begin(), 2)), diploidOpt.failedSampleFTLabel);
    
    // Designed the case-6 and case-7. MinGTScore = 50 
    diploidOpt.minPassGTScore = 50;
    diploidInfo.filters.clear();
    scoreDiploidSV(diploidOpt, scanner.operator*(), diploidDeriv, dFilter, junctionData, diploidInfo);
    BOOST_REQUIRE_EQUAL(diploidInfo.samples[0].filters.size(), 2);
    BOOST_REQUIRE_EQUAL(*(diploidInfo.samples[0].filters.begin()), diploidOpt.homRefLabel);
    BOOST_REQUIRE_EQUAL(*(std::next(diploidInfo.samples[0].filters.begin(), 1)), diploidOpt.minGTFilterLabel);
    
    // Designed the case-5
    const std::shared_ptr<TestChromosomeDepthFileMaker> depthFileMaker(new TestChromosomeDepthFileMaker());
    const std::string depthFileName(depthFileMaker.operator*().getFilename());
    buildTestChromosomeDepthFile(depthFileName);
    ChromDepthFilterUtil depthFilterUtil2(depthFileName, 2, bamHeader);
    SVScoreInfo scoreInfo2;
    scoreInfo2.bp2MaxDepth = 120;
    scoreInfo2.samples.resize(1);
    scoreInfo2.samples.push_back(diploidSampleInfo);
    junctionData.clear();
    diploidInfo.filters.clear();
    callInfo.init(candidate, evidence, scoreInfo2, 0.2);
    junctionData.push_back(callInfo);
    scoreDiploidSV(diploidOpt, scanner.operator*(), diploidDeriv, depthFilterUtil2, junctionData, diploidInfo);
    BOOST_REQUIRE_EQUAL(diploidInfo.filters.size(), 4);
    BOOST_REQUIRE_EQUAL(*(diploidInfo.filters.begin()), diploidOpt.maxDepthFilterLabel);
    BOOST_REQUIRE_EQUAL(*(std::next(diploidInfo.filters.begin(), 1)), diploidOpt.minAltFilterLabel);
    BOOST_REQUIRE_EQUAL(*(std::next(diploidInfo.filters.begin(), 2)), diploidOpt.noPairSupportLabel);
    BOOST_REQUIRE_EQUAL(*(std::next(diploidInfo.filters.begin(), 3)), diploidOpt.failedSampleFTLabel);
}

// Test the following cases:
// 1. when MQ0 fraction of BP1 or BP2 is greater than control filtration based on MQ0 fraction (0.4)
// 2. when max depth of BP1 or BP2 is greater than chromosome depth
BOOST_AUTO_TEST_CASE( test_TumorSVLabel )
{
    const bam_header_info bamHeader(buildTestBamHeader());
    CallOptionsTumor optionsTumor;
    // Dummy depth file that means depth file is not present
    ChromDepthFilterUtil depthFilterUtil1("", 0.2, bamHeader);
    SVCandidate candidate;
    candidate.bp1.interval = GenomeInterval(0, 100, 200);
    candidate.bp2.interval = GenomeInterval(0, 400, 500);
    candidate.insertSeq = "AGCTGACTGATCAGT";
    SVScoreInfoTumor scoreInfoTumor;
    SVScoreInfo scoreInfo1;
    scoreInfo1.bp1MQ0Frac = 0.2f;
    scoreInfo1.bp2MQ0Frac = 0.3f;
    SVEvidence evidence;
    JunctionCallInfo callInfo;
    callInfo.init(candidate, evidence, scoreInfo1, 0.2);
    std::vector <JunctionCallInfo> callInfos;
    callInfos.push_back(callInfo);
    scoreTumorSV(optionsTumor, depthFilterUtil1, callInfos, scoreInfoTumor);
    // All the cases are fine. No label will be added.
    BOOST_REQUIRE_EQUAL(scoreInfoTumor.filters.size(), 0);
    
    // Designed the case-1 for BP1
    callInfos.clear();
    scoreInfo1.bp1MQ0Frac = 0.2f;
    scoreInfo1.bp2MQ0Frac = 0.5f;
    callInfo.init(candidate, evidence, scoreInfo1, 0.2);
    callInfos.push_back(callInfo);
    scoreTumorSV(optionsTumor, depthFilterUtil1, callInfos, scoreInfoTumor);
    BOOST_REQUIRE_EQUAL(scoreInfoTumor.filters.size(), 1);
    BOOST_REQUIRE_EQUAL(*(scoreInfoTumor.filters.begin()), optionsTumor.maxMQ0FracLabel);
    
    // Designed the case-1 for BP2
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
    
    // Designed the case-2 
    callInfos.clear();
    scoreInfoTumor.filters.clear();
    const std::shared_ptr<TestChromosomeDepthFileMaker> depthFileMaker(new TestChromosomeDepthFileMaker());
    const std::string depthFileName(depthFileMaker.operator*().getFilename());
    buildTestChromosomeDepthFile(depthFileName);
    ChromDepthFilterUtil depthFilterUtil2(depthFileName, 2, bamHeader);
    SVScoreInfo scoreInfo2;
    scoreInfo2.bp2MaxDepth = 120;
    callInfo.init(candidate, evidence, scoreInfo2, 0.2);
    callInfos.push_back(callInfo);
    scoreTumorSV(optionsTumor, depthFilterUtil2, callInfos, scoreInfoTumor);
    BOOST_REQUIRE_EQUAL(scoreInfoTumor.filters.size(), 1);
    BOOST_REQUIRE_EQUAL(*(scoreInfoTumor.filters.begin()), optionsTumor.maxDepthFilterLabel);
}

// Following cases need to be tested:
// 1. When SV is imprecise. Add Imrecise label.
// 2. When difference between two breakpoints is less than Min length for passing fusions (100000).
//    Add a local label.
// 3. when split read count or confident spanning pair count is 0. Add a rnaFilterLabel.
BOOST_AUTO_TEST_CASE( test_scoreRNASV )
{
    SVSampleInfo sampleInfo;
    SVCandidate candidate;
    SVScoreInfoRna scoreInfoRna;
    scoreRNASV(sampleInfo, candidate, scoreInfoRna);
    // Designed the case-1
    BOOST_REQUIRE_EQUAL(scoreInfoRna.filters.size(), 1);
    BOOST_REQUIRE_EQUAL(*(scoreInfoRna.filters.begin()), SVScoreInfoRna::impreciseLabel);
    
    // Designed the case-2
    candidate.setPrecise();
    scoreInfoRna.filters.clear();
    candidate.bp1.interval = GenomeInterval(0, 100, 200);
    candidate.bp2.interval = GenomeInterval(0, 400, 450);
    sampleInfo.alt.splitReadCount = 1;
    sampleInfo.alt.confidentSpanningPairCount = 1;
    scoreRNASV(sampleInfo, candidate, scoreInfoRna);
    BOOST_REQUIRE_EQUAL(scoreInfoRna.filters.size(), 1);
    BOOST_REQUIRE_EQUAL(*(scoreInfoRna.filters.begin()), SVScoreInfoRna::localLabel);
    
    // Designed the case-3 when splitReadCount = 0 and confidentSpanningPairCount = 1
    scoreInfoRna.filters.clear();
    candidate.bp1.interval = GenomeInterval(0, 100, 200);
    candidate.bp2.interval = GenomeInterval(1, 400, 450);
    sampleInfo.alt.splitReadCount = 0;
    sampleInfo.alt.confidentSpanningPairCount = 1;
    scoreRNASV(sampleInfo, candidate, scoreInfoRna);
    BOOST_REQUIRE_EQUAL(scoreInfoRna.filters.size(), 1);
    BOOST_REQUIRE_EQUAL(*(scoreInfoRna.filters.begin()), SVScoreInfoRna::rnaFilterLabel);
    
    // Designed the case-3 when splitReadCount = 1 and confidentSpanningPairCount = 0
    scoreInfoRna.filters.clear();
    candidate.bp1.interval = GenomeInterval(0, 100, 200);
    candidate.bp2.interval = GenomeInterval(1, 400, 450);
    sampleInfo.alt.splitReadCount = 1;
    sampleInfo.alt.confidentSpanningPairCount = 0;
    scoreRNASV(sampleInfo, candidate, scoreInfoRna);
    BOOST_REQUIRE_EQUAL(scoreInfoRna.filters.size(), 1);
    BOOST_REQUIRE_EQUAL(*(scoreInfoRna.filters.begin()), SVScoreInfoRna::rnaFilterLabel);
    
    // Designed the case-3 when splitReadCount = 0 and confidentSpanningPairCount = 0
    scoreInfoRna.filters.clear();
    candidate.bp1.interval = GenomeInterval(0, 100, 200);
    candidate.bp2.interval = GenomeInterval(0, 400, 450);
    sampleInfo.alt.splitReadCount = 0;
    sampleInfo.alt.confidentSpanningPairCount = 0;
    scoreRNASV(sampleInfo, candidate, scoreInfoRna);
    BOOST_REQUIRE_EQUAL(scoreInfoRna.filters.size(), 2);
    BOOST_REQUIRE_EQUAL(*(scoreInfoRna.filters.begin()), SVScoreInfoRna::localLabel);
    BOOST_REQUIRE_EQUAL(*(std::next(scoreInfoRna.filters.begin(), 1)), SVScoreInfoRna::rnaFilterLabel);
}

// If tier2 support
//     spanning pair count = confidentSemiMappedSpanningPairCount * spanning weight
// else
//     spanning pair count = confidentSpanningPairCount * spanning weight
BOOST_AUTO_TEST_CASE( test_getSpanningPairCount )
{
    SVSampleAlleleInfo alleleInfo;
    alleleInfo.confidentSpanningPairCount = 100;
    alleleInfo.confidentSemiMappedSpanningPairCount = 50;
    BOOST_REQUIRE_EQUAL(getSpanningPairCount(alleleInfo, 0.2, false), 20);
    // for tier2 support
    BOOST_REQUIRE_EQUAL(getSpanningPairCount(alleleInfo, 0.2, true), 10);
}

// If tier2 support
//     spanning pair count = confidentSplitReadCount + confidentSemiMappedSpanningPairCount * spanning weight
// else
//     spanning pair count = confidentSplitReadCount + confidentSpanningPairCount * spanning weight
BOOST_AUTO_TEST_CASE( test_getSupportCount )
{
    SVSampleAlleleInfo alleleInfo;
    alleleInfo.confidentSpanningPairCount = 100;
    alleleInfo.confidentSemiMappedSpanningPairCount = 50;
    alleleInfo.confidentSplitReadCount = 10;
    BOOST_REQUIRE_EQUAL(getSupportCount(alleleInfo, 0.2, false), 30);
    BOOST_REQUIRE_EQUAL(getSupportCount(alleleInfo, 0.2, true), 20);
}

// Test the fraction of allele support count to total support
// counts (allle support + ref support)
BOOST_AUTO_TEST_CASE( test_estimateSomaticMutationFreq )
{
    SVCandidate candidate;
    SVEvidence evidence;
    
    // All counts are zero so fraction should be also 0.
    std::vector<JunctionCallInfo> callInfos;
    callInfos.resize(1);
    SVScoreInfo scoreInfo1;
    SVSampleInfo sampleInfo;
    scoreInfo1.samples.resize(1);
    scoreInfo1.samples.push_back(sampleInfo);
    callInfos[0].init(candidate, evidence, scoreInfo1, 0.2);
    BOOST_REQUIRE_EQUAL(estimateSomaticMutationFreq(0, callInfos, false), 0);
    
    // allele support count > 0 and ref support count is 0. So fraction of 
    // allele support = 1
    sampleInfo.alt.confidentSpanningPairCount = 100;
    sampleInfo.alt.confidentSemiMappedSpanningPairCount = 50;
    sampleInfo.alt.confidentSplitReadCount = 10;
    SVScoreInfo scoreInfo2;
    scoreInfo2.samples.push_back(sampleInfo);
    callInfos[0].init(candidate, evidence, scoreInfo2, 0.2);
    BOOST_REQUIRE_EQUAL(estimateSomaticMutationFreq(0, callInfos, false), 1);

    // allele support count = 0 and ref support count > 0. So fraction of 
    // allele support = 0
    sampleInfo.alt.confidentSpanningPairCount = 0;
    sampleInfo.alt.confidentSemiMappedSpanningPairCount = 0;
    sampleInfo.alt.confidentSplitReadCount = 0;
    sampleInfo.ref.confidentSpanningPairCount = 100;
    sampleInfo.ref.confidentSemiMappedSpanningPairCount = 50;
    sampleInfo.ref.confidentSplitReadCount = 10;
    SVScoreInfo scoreInfo3;
    scoreInfo3.samples.push_back(sampleInfo);
    callInfos[0].init(candidate, evidence, scoreInfo3, 0.2);
    BOOST_REQUIRE_EQUAL(estimateSomaticMutationFreq(0, callInfos, false), 0);

    // allele support count > 0 and ref support count is >. So fraction of 
    // allele support = allele support / (allele support + ref support)
    sampleInfo.alt.confidentSpanningPairCount = 100;
    sampleInfo.alt.confidentSemiMappedSpanningPairCount = 50;
    sampleInfo.alt.confidentSplitReadCount = 10;
    sampleInfo.ref.confidentSpanningPairCount = 100;
    sampleInfo.ref.confidentSemiMappedSpanningPairCount = 50;
    sampleInfo.ref.confidentSplitReadCount = 10;
    SVScoreInfo scoreInfo4;
    scoreInfo4.samples.push_back(sampleInfo);
    callInfos[0].init(candidate, evidence, scoreInfo4, 0.2);
    BOOST_REQUIRE_EQUAL(estimateSomaticMutationFreq(0, callInfos, false), 0.5);
}

// This is applicable for somatic caller.
// Test the fraction of allele support count (normal + tumor) to total support
// counts (allle support + ref support of normal and tumor)
BOOST_AUTO_TEST_CASE( test_estimateNoiseMutationFreq )
{
    SVCandidate candidate;
    SVEvidence evidence;
    // All counts are zero for both normal and tumor
    std::vector<JunctionCallInfo> callInfos;
    callInfos.resize(1);
    SVScoreInfo scoreInfo1;
    SVSampleInfo normalSampleInfo;
    SVSampleInfo tumorSampleInfo;
    scoreInfo1.samples.resize(2);
    scoreInfo1.samples.push_back(normalSampleInfo);
    scoreInfo1.samples.push_back(tumorSampleInfo);
    callInfos[0].init(candidate, evidence, scoreInfo1, 0.2);
    BOOST_REQUIRE_EQUAL(estimateNoiseMutationFreq(0, 1, callInfos, false), 0);
    
    // allele support count (tumor + normal) > 0 and ref support count is 0. So fraction of 
    // allele support = 1
    normalSampleInfo.alt.confidentSpanningPairCount = 100;
    normalSampleInfo.alt.confidentSemiMappedSpanningPairCount = 50;
    normalSampleInfo.alt.confidentSplitReadCount = 10;
    tumorSampleInfo.alt.confidentSpanningPairCount = 100;
    tumorSampleInfo.alt.confidentSemiMappedSpanningPairCount = 50;
    tumorSampleInfo.alt.confidentSplitReadCount = 10;
    SVScoreInfo scoreInfo2;
    scoreInfo2.samples.push_back(normalSampleInfo);
    scoreInfo2.samples.push_back(tumorSampleInfo);
    callInfos[0].init(candidate, evidence, scoreInfo2, 0.2);
    BOOST_REQUIRE_EQUAL(estimateNoiseMutationFreq(0, 1, callInfos, false), 1);

    // allele support count (tumor + normal) = 0 and ref support count
    // (tumor + normal) > 0. So fraction of allele support = 0
    normalSampleInfo.alt.confidentSpanningPairCount = 0;
    normalSampleInfo.alt.confidentSemiMappedSpanningPairCount = 0;
    normalSampleInfo.alt.confidentSplitReadCount = 0;
    normalSampleInfo.ref.confidentSpanningPairCount = 100;
    normalSampleInfo.ref.confidentSemiMappedSpanningPairCount = 50;
    normalSampleInfo.ref.confidentSplitReadCount = 10;
    tumorSampleInfo.alt.confidentSpanningPairCount = 0;
    tumorSampleInfo.alt.confidentSemiMappedSpanningPairCount = 0;
    tumorSampleInfo.alt.confidentSplitReadCount = 0;
    tumorSampleInfo.ref.confidentSpanningPairCount = 100;
    tumorSampleInfo.ref.confidentSemiMappedSpanningPairCount = 50;
    tumorSampleInfo.ref.confidentSplitReadCount = 10;
    SVScoreInfo scoreInfo3;
    scoreInfo3.samples.push_back(normalSampleInfo);
    scoreInfo3.samples.push_back(tumorSampleInfo);
    callInfos[0].init(candidate, evidence, scoreInfo3, 0.2);
    BOOST_REQUIRE_EQUAL(estimateNoiseMutationFreq(0, 1, callInfos, false), 0);

    // allele support count (tumor + normal) > 0 and ref support count (tumor + normal) > 0.
    // So fraction of allele support = allele support / (allele support + ref support)
    normalSampleInfo.alt.confidentSpanningPairCount = 100;
    normalSampleInfo.alt.confidentSemiMappedSpanningPairCount = 50;
    normalSampleInfo.alt.confidentSplitReadCount = 10;
    normalSampleInfo.ref.confidentSpanningPairCount = 100;
    normalSampleInfo.ref.confidentSemiMappedSpanningPairCount = 50;
    normalSampleInfo.ref.confidentSplitReadCount = 10;
    tumorSampleInfo.alt.confidentSpanningPairCount = 100;
    tumorSampleInfo.alt.confidentSemiMappedSpanningPairCount = 50;
    tumorSampleInfo.alt.confidentSplitReadCount = 10;
    tumorSampleInfo.ref.confidentSpanningPairCount = 100;
    tumorSampleInfo.ref.confidentSemiMappedSpanningPairCount = 50;
    tumorSampleInfo.ref.confidentSplitReadCount = 10;
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
// So if refLnFragLhood of ref allele is R and altLnFragLhood of alt allele
// is L, then score(S) of a genotype (gt = REF or HET or HOM or n or s) is
// S = log_sum(R + lnPreComp(gt), L + LnPre(gt)), where,
//     LnPre(gt) = log(Prior probability of expected alt allele)
//     lnPreComp(gt) = log(1-Prior probability of expected alt allele)
//     log_sum(x1,x2) = x1 + log(1+exp(x2-x1)) if x2>x1
//                    = x2 + log(1+exp(x1-x2)) if x1>x2
BOOST_AUTO_TEST_CASE( test_computeSomaticSampleLoghood )
{
    float spanningPairWeight = 0.2;
    ProbSet refChimeraProb(0.2);
    ProbSet altChimeraProb(0.3);
    ProbSet refSplitMapProb(0.1);
    ProbSet altSplitMapProb(0.5);
    SVEvidence::evidenceTrack_t evidenceTrack;
    std::string fragLabel = "frag-1";
    SVFragmentEvidence fragmentEvidence1;
    fragmentEvidence1.alt.bp1.isFragmentSupport = true;
    fragmentEvidence1.alt.bp1.fragLengthProb = 0.4;
    fragmentEvidence1.ref.bp1.read1.isSplitSupport = true;
    fragmentEvidence1.ref.bp1.read1.splitLnLhood = 0.2;
    evidenceTrack[fragLabel] = fragmentEvidence1;
    std::array<double,SOMATIC_GT::SIZE> loglhood;
    for (unsigned gt(0); gt < SOMATIC_GT::SIZE; ++gt)
    {
        loglhood[gt] = 0;
    }
    // fragment was not evaluated for pair or split support for either allele
    computeSomaticSampleLoghood(spanningPairWeight, evidenceTrack, 0.4, 0.25, false, false, refChimeraProb, altChimeraProb, refSplitMapProb, altSplitMapProb, loglhood);
    for (unsigned gt(0); gt < SOMATIC_GT::SIZE; ++gt)
    {
        BOOST_REQUIRE_EQUAL(loglhood[gt], 0);
    }
    
    // fragment was evaluated for pair or split support for either allele
    fragmentEvidence1.read1.isScanned = true;
    fragmentEvidence1.read2.isScanned = true;
    fragmentEvidence1.read1.setAnchored(true);
    fragmentEvidence1.read2.setAnchored(true);
    fragmentEvidence1.alt.bp1.isFragmentSupport = true;
    evidenceTrack[fragLabel] = fragmentEvidence1;
    computeSomaticSampleLoghood(spanningPairWeight, evidenceTrack, 0.4, 0.25, false, false, refChimeraProb, altChimeraProb, refSplitMapProb, altSplitMapProb, loglhood);
    // above mentioned formula is applied in loglhood
    for (unsigned gt(0); gt < SOMATIC_GT::SIZE; ++gt)
    {
        BOOST_REQUIRE(loglhood[gt] != 0); // just checking of non-zero because of floating point issue
    }

}

// Following cases need to be tested
// 1. When somatic variant score is less than minimum somatic quality(30). Add minScoreLabel.
// 2. when max depth of BP1 or BP2 is greater than chromosome depth. Add maxDepth Label.
// 3. when MQ0 fraction of BP1 or BP2 is greater than control filtration based on MQ0 fraction (0.4). Add
//    maxMQ0FracLabel
BOOST_AUTO_TEST_CASE( test_SomaticSVLabel )
{
    // somatic caller options
    CallOptionsSomatic somaticOpt;
    CallOptionsSomaticDeriv somaticDopt(somaticOpt);
    // created chromsome depth file
    const std::shared_ptr<TestChromosomeDepthFileMaker> depthFileMaker(new TestChromosomeDepthFileMaker());
    const std::string depthFileName(depthFileMaker.operator*().getFilename());
    buildTestChromosomeDepthFile(depthFileName);
    const bam_header_info bamHeader(buildTestBamHeader());
    ChromDepthFilterUtil depthFilterUtil(depthFileName, 2, bamHeader);

    std::vector<JunctionCallInfo> junctionData;
    junctionData.resize(1); // adding one junction
    // Adding few counts to tumor and normal sample
    SVSampleInfo normalSampleInfo;
    SVSampleInfo tumorSampleInfo;
    normalSampleInfo.alt.confidentSpanningPairCount = 100;
    normalSampleInfo.alt.confidentSemiMappedSpanningPairCount = 50;
    normalSampleInfo.alt.confidentSplitReadCount = 10;
    normalSampleInfo.ref.confidentSpanningPairCount = 100;
    normalSampleInfo.ref.confidentSemiMappedSpanningPairCount = 50;
    normalSampleInfo.ref.confidentSplitReadCount = 10;
    tumorSampleInfo.alt.confidentSpanningPairCount = 100;
    tumorSampleInfo.alt.confidentSemiMappedSpanningPairCount = 50;
    tumorSampleInfo.alt.confidentSplitReadCount = 10;
    tumorSampleInfo.ref.confidentSpanningPairCount = 100;
    tumorSampleInfo.ref.confidentSemiMappedSpanningPairCount = 50;
    tumorSampleInfo.ref.confidentSplitReadCount = 10;
    // somatic SV score
    SVScoreInfo scoreInfo;
    scoreInfo.samples.push_back(normalSampleInfo);
    scoreInfo.samples.push_back(tumorSampleInfo);
    SVEvidence evidence;
    evidence.samples.resize(2);
    SVCandidate candidate;
    junctionData[0].init(candidate, evidence, scoreInfo, 0.2);
    SVScoreInfoSomatic somaticInfo;
    somaticInfo.somaticScore = 20;
    // Designed the case-1
    scoreSomaticSV(2, 1, somaticOpt, somaticDopt, depthFilterUtil, junctionData, somaticInfo);
    BOOST_REQUIRE_EQUAL(somaticInfo.filters.size(), 1);
    BOOST_REQUIRE_EQUAL(*(somaticInfo.filters.begin()), somaticOpt.minSomaticScoreLabel);

    // Designed the case-1 and case-2
    scoreInfo.bp2MaxDepth = 120;
    somaticInfo.filters.clear();
    somaticInfo.somaticScore = 40;
    junctionData[0].init(candidate, evidence, scoreInfo, 0.2);
    scoreSomaticSV(2, 1, somaticOpt, somaticDopt, depthFilterUtil, junctionData, somaticInfo);
    BOOST_REQUIRE_EQUAL(somaticInfo.filters.size(), 2);
    BOOST_REQUIRE_EQUAL(*(somaticInfo.filters.begin()), somaticOpt.maxDepthFilterLabel);
    BOOST_REQUIRE_EQUAL(*(std::next(somaticInfo.filters.begin(), 1)), somaticOpt.minSomaticScoreLabel);

    // Designed the case-1 and case-3
    scoreInfo.bp2MaxDepth = 20;
    scoreInfo.bp2MQ0Frac = 0.6;
    somaticInfo.filters.clear();
    somaticInfo.somaticScore = 40;
    junctionData[0].init(candidate, evidence, scoreInfo, 0.2);
    scoreSomaticSV(2, 1, somaticOpt, somaticDopt, depthFilterUtil, junctionData, somaticInfo);
    BOOST_REQUIRE_EQUAL(somaticInfo.filters.size(), 2);
    BOOST_REQUIRE_EQUAL(*(somaticInfo.filters.begin()), somaticOpt.maxMQ0FracLabel);
    BOOST_REQUIRE_EQUAL(*(std::next(somaticInfo.filters.begin(), 1)), somaticOpt.minSomaticScoreLabel);
}

// Test the fraction of zero mapping quality reads
BOOST_AUTO_TEST_CASE( test_getBreakendMaxMappedDepthAndMQ0 )
{
    const bam_header_info bamHeader(buildTestBamHeader());
    std::unique_ptr<SVLocusScanner> scanner(buildTestSVLocusScanner(bamHeader));
    GSCOptions options;
    options.alignFileOpt.alignmentFilenames = {bamFilename};
    SVScorer scorer(options, scanner.operator*(), bamHeader);
    TestSVScorer fSVScorer;
    SVBreakend breakend;
    breakend.interval = GenomeInterval(0, 109, 110);
    breakend.state = SVBreakendState::RIGHT_OPEN;
    float mq0Frac = 0;
    unsigned maxDepth(1000);
    // No depth cutoff in this test cases 
    // Total 18 reads have mapping quality 0 out of 20 reads. So fraction is 18/20 = 0.9
    fSVScorer.getBreakendMaxMappedDepthAndMQ0(scorer, false, false, 100, breakend, maxDepth, mq0Frac);
    BOOST_REQUIRE_EQUAL(((int)(mq0Frac * 100 + 0.5)) / 100.0, 0.9);
    
    // depth cutoff is 12. If depth of a location more than 12, it will return the fraction.
    // So here total 11 reads have mapping quality 0 out of 13 reads. So fraction is 11/13 = ~0.85.
    fSVScorer.getBreakendMaxMappedDepthAndMQ0(scorer, false, true, 12, breakend, maxDepth, mq0Frac);
    BOOST_REQUIRE_EQUAL(((int)(mq0Frac * 100 + 0.5)) / 100.0, 0.85);
}

// Test Whether a specific scoring model is performed based on arguments
BOOST_AUTO_TEST_CASE( test_computeAllScoreModel )
{
    const bam_header_info bamHeader(buildTestBamHeader());
    std::unique_ptr<SVLocusScanner> scanner(buildTestSVLocusScanner(bamHeader));
    GSCOptions options;
    options.alignFileOpt.alignmentFilenames = {bamFilename};
    SVScorer scorer(options, scanner.operator*(), bamHeader);
    TestSVScorer fSVScorer;
    SVModelScoreInfo modelScoreInfo;
    modelScoreInfo.setSampleCount(1, 1);

    SVCandidate candidate;
    candidate.bp1.interval = GenomeInterval(0, 100, 200);
    candidate.bp2.interval = GenomeInterval(0, 400, 500);
    candidate.insertSeq = "AGCTGACTGATCAGT";
    SVScoreInfo scoreInfo;
    scoreInfo.bp1MQ0Frac = 0.2f;
    scoreInfo.bp2MQ0Frac = 0.5f;
    SVEvidence evidence;
    JunctionCallInfo callInfo;
    callInfo.init(candidate, evidence, scoreInfo, 0.2);
    std::vector <JunctionCallInfo> callInfos;
    callInfos.push_back(callInfo);
    // Tumor scoring model
    fSVScorer.computeAllScoreModels(scorer, false, true, callInfos, modelScoreInfo);
    BOOST_REQUIRE_EQUAL(modelScoreInfo.tumor.filters.size(), 1);
    callInfos.clear();
    candidate.setPrecise();
    candidate.bp1.interval = GenomeInterval(0, 100, 200);
    candidate.bp2.interval = GenomeInterval(0, 400, 450);
    modelScoreInfo.base.samples[0].alt.splitReadCount = 1;
    modelScoreInfo.base.samples[0].alt.confidentSpanningPairCount = 1;
    SVFragmentEvidence fragmentEvidence;
    fragmentEvidence.read1.isScanned = true;
    fragmentEvidence.read2.isScanned = true;
    fragmentEvidence.read1.setAnchored(true);
    fragmentEvidence.read2.setAnchored(true);
    fragmentEvidence.alt.bp1.isFragmentSupport = true;
    SVEvidence::evidenceTrack_t evidenceTrack;
    std::string fragLabel = "frag-1";
    evidenceTrack[fragLabel] = fragmentEvidence;
    evidence.samples.push_back(evidenceTrack);
    scoreInfo.samples.resize(1);
    callInfo.init(candidate, evidence, scoreInfo, 0.2);
    callInfos.push_back(callInfo);
    
    // Only Diploid scoring model
    fSVScorer.computeAllScoreModels(scorer, false, false, callInfos, modelScoreInfo);
    BOOST_REQUIRE_EQUAL(modelScoreInfo.diploid.filters.size(), 4);

    // somatic and diploid scoring model
    CallOptionsSomatic somaticOpt;
    CallOptionsSomaticDeriv somaticDopt(somaticOpt);
    const std::shared_ptr<TestChromosomeDepthFileMaker> depthFileMaker(new TestChromosomeDepthFileMaker());
    const std::string depthFileName(depthFileMaker.operator*().getFilename());
    buildTestChromosomeDepthFile(depthFileName);
    ChromDepthFilterUtil depthFilterUtil(depthFileName, 2, bamHeader);
    std::vector<JunctionCallInfo> junctionData;
    junctionData.resize(1);
    SVSampleInfo normalSampleInfo;
    SVSampleInfo tumorSampleInfo;
    normalSampleInfo.alt.confidentSpanningPairCount = 100;
    normalSampleInfo.alt.confidentSemiMappedSpanningPairCount = 50;
    normalSampleInfo.alt.confidentSplitReadCount = 10;
    normalSampleInfo.ref.confidentSpanningPairCount = 100;
    normalSampleInfo.ref.confidentSemiMappedSpanningPairCount = 50;
    normalSampleInfo.ref.confidentSplitReadCount = 10;
    tumorSampleInfo.alt.confidentSpanningPairCount = 100;
    tumorSampleInfo.alt.confidentSemiMappedSpanningPairCount = 50;
    tumorSampleInfo.alt.confidentSplitReadCount = 10;
    tumorSampleInfo.ref.confidentSpanningPairCount = 100;
    tumorSampleInfo.ref.confidentSemiMappedSpanningPairCount = 50;
    tumorSampleInfo.ref.confidentSplitReadCount = 10;
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
    fSVScorer.computeAllScoreModels(scorer, true, false, junctionData, modelScoreInfo);
    BOOST_REQUIRE_EQUAL(modelScoreInfo.somatic.filters.size(), 1);

    // RNA scoring model
    modelScoreInfo.setSampleCount(1, 1);
    GSCOptions options2;
    options2.alignFileOpt.alignmentFilenames = {bamFilename};
    options2.isRNA = true;
    SVScorer scorer2(options2, scanner.operator*(), bamHeader);
    TestSVScorer fSVScorer2;
    fSVScorer2.computeAllScoreModels(scorer2, false, false, callInfos, modelScoreInfo);
    BOOST_REQUIRE_EQUAL(modelScoreInfo.rna.filters.size(), 2);
}

// Few minor logic checking have been done in this test case
// 1. When no valid junction is present. It will thow an exception
// 2. When only one unfiltered junction is present
// 3. When more than one unfiltered junction is present
BOOST_AUTO_TEST_CASE( test_ScoreSV )
{
    // Designed case-1
    SVCandidateSetData svData;
    std::vector<SVCandidateAssemblyData> mjAssemblyData;
    SVMultiJunctionCandidate mjSV;
    std::vector<SVId> mjSVId;
    std::vector<bool> isJunctionFiltered;
    std::vector<SVModelScoreInfo> mjModelScoreInfo;
    SVModelScoreInfo mjJointModelScoreInfo;
    bool isMJEvent;
    SupportSamples svSupports;
    const bam_header_info bamHeader(buildTestBamHeader());
    std::unique_ptr<SVLocusScanner> scanner(buildTestSVLocusScanner(bamHeader));
    GSCOptions options;
    options.alignFileOpt.alignmentFilenames = {bamFilename};
    SVScorer scorer(options, scanner.operator*(), bamHeader);
    BOOST_CHECK_THROW(scorer.scoreSV(svData, mjAssemblyData, mjSV, mjSVId, isJunctionFiltered,
                                     false, true, mjModelScoreInfo, mjJointModelScoreInfo, isMJEvent, svSupports),
                                     illumina::common::GeneralException);
    
    // Designed case-2. Here number of unfiltered multi-junction event is 1.
    // So isMJEvent should be false here.
    SVCandidate candidate1;
    candidate1.bp1.interval = GenomeInterval(0, 100, 200);
    candidate1.bp2.interval = GenomeInterval(0, 400, 500);
    candidate1.bp1.state = SVBreakendState::RIGHT_OPEN;
    candidate1.bp2.state = SVBreakendState::LEFT_OPEN;
    candidate1.insertSeq = "AGCTGACTGATCAGT";
    candidate1.setPrecise();
    candidate1.assemblyAlignIndex = 0;
    mjSV.junction.push_back(candidate1);
    SVId id1;
    SVId id2;
    mjSVId.push_back(id1);
    mjSVId.push_back(id2);
    isJunctionFiltered.push_back(false);
    SVModelScoreInfo modelScoreInfo;
    modelScoreInfo.setSampleCount(1, 1);
    mjModelScoreInfo.push_back(modelScoreInfo);
    mjJointModelScoreInfo.setSampleCount(1, 1);
    TestSVScorer fSVScorer;
    fSVScorer.setSampleCount(scorer, 1, 1);
    svSupports.supportSamples.resize(1);

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
    SVCandidateSetSequenceFragmentSampleGroup& group = svData.getDataGroup(0);
    group.add(bamHeader, bamRecord1, false, true, true);
    group.add(bamHeader, bamRecord2, false, true, true);
    SVSequenceFragmentAssociation association(0, SVEvidenceType::PAIR);
    group.begin().operator*().svLink.push_back(association);

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
    mjAssemblyData.push_back(candidateAssemblyData);
    scorer.scoreSV(svData, mjAssemblyData, mjSV, mjSVId, isJunctionFiltered,
                                     false, true, mjModelScoreInfo, mjJointModelScoreInfo, isMJEvent, svSupports);
    BOOST_REQUIRE(!isMJEvent);
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
