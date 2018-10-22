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
#include "manta/SVLocusEvidenceCount.hh"
#include "test/testAlignmentDataUtil.hh"
#include "test/testSVLocusScanner.hh"
#include "test/testUtil.hh"
#include "test/testFileMakers.hh"

#include "SVFinder.hh"
#include "SVFinder.cpp"

BOOST_AUTO_TEST_SUITE( SVFinderTest_test_suite )


BOOST_AUTO_TEST_CASE( test_SpanningNoiseRate )
{
    AllSampleReadCounts counts;
    counts.setSampleCount(2);
    SampleReadCounts sample1(counts.getSampleCounts(0));
    sample1.input.evidenceCount.anom = 10;
    sample1.input.evidenceCount.split = 5;
    sample1.input.evidenceCount.anomAndSplit = 4;
    sample1.input.evidenceCount.total = 19;
    counts.getSampleCounts(0).merge(sample1);

    SampleReadCounts sample2(counts.getSampleCounts(1));
    sample2.input.evidenceCount.anom = 25;
    sample2.input.evidenceCount.split = 5;
    sample2.input.evidenceCount.anomAndSplit = 10;
    sample2.input.evidenceCount.total = 40;
    counts.getSampleCounts(1).merge(sample2);

    BOOST_REQUIRE_EQUAL(getSpanningNoiseRate(counts, 0), 0.020608439646712464);
    BOOST_REQUIRE_EQUAL(getSpanningNoiseRate(counts, 1), 0.028846153846153848);
}

BOOST_AUTO_TEST_CASE( test_AssemblyNoiseRate )
{
    AllSampleReadCounts counts;
    counts.setSampleCount(2);
    SampleReadCounts sample1(counts.getSampleCounts(0));
    sample1.input.evidenceCount.assm = 10;
    sample1.input.evidenceCount.total = 19;
    counts.getSampleCounts(0).merge(sample1);

    SampleReadCounts sample2(counts.getSampleCounts(1));
    sample2.input.evidenceCount.assm = 25;
    sample2.input.evidenceCount.total = 40;
    counts.getSampleCounts(1).merge(sample2);

    BOOST_REQUIRE_EQUAL(getAssemblyNoiseRate(counts, 0), 0.019627085377821395);
    BOOST_REQUIRE_EQUAL(getAssemblyNoiseRate(counts, 1), 0.033653846153846152);
}

BOOST_AUTO_TEST_CASE( test_AddSVNode )
{
    // test if read supports an SV on this edge, if so, add to SVData

    const bam_header_info bamHeader(buildTestBamHeader());
    std::unique_ptr<SVLocusScanner> scanner(buildTestSVLocusScanner(bamHeader));
    SampleEvidenceCounts eCounts;

    const unsigned defaultReadGroupIndex = 0;
    reference_contig_segment seq = reference_contig_segment();
    const reference_contig_segment& refSeq(seq);

    bam_record supplementSASplitRead;
    buildTestBamRecord(supplementSASplitRead);
    addSupplementaryAlignmentEvidence(supplementSASplitRead);

    // large insertion is SV evidence
    bam_record largeInsertionRead;
    buildTestBamRecord(largeInsertionRead, 0, 200, 0, 300, 100, 15, "100M2000I100M");
    largeInsertionRead.set_qname("large_insertion");

    SVLocus locus1;
    locus1.addNode(GenomeInterval(0,80,120));
    locus1.addNode(GenomeInterval(0,279,319));
    locus1.addNode(GenomeInterval(0,410,450));

    SVCandidateSetSequenceFragmentSampleGroup svDatagroup;

    // read overlaps
    addSVNodeRead(bamHeader, scanner.operator*(), static_cast<const SVLocus&>(locus1).getNode(0), static_cast<const SVLocus&>(locus1).getNode(1),
            supplementSASplitRead, defaultReadGroupIndex, true, refSeq, true, false, svDatagroup, eCounts);

    BOOST_REQUIRE_EQUAL(svDatagroup.size(), 1u);

    // read does not overlap
    BOOST_REQUIRE_EQUAL(svDatagroup.size(), 1u);
    addSVNodeRead(bamHeader, scanner.operator*(), static_cast<const SVLocus&>(locus1).getNode(0), static_cast<const SVLocus&>(locus1).getNode(1),
            largeInsertionRead, defaultReadGroupIndex, true, refSeq, true, false, svDatagroup, eCounts);
    // read overlaps
    addSVNodeRead(bamHeader, scanner.operator*(), static_cast<const SVLocus&>(locus1).getNode(2), static_cast<const SVLocus&>(locus1).getNode(1),
            largeInsertionRead, defaultReadGroupIndex, true, refSeq, true, false, svDatagroup, eCounts);
    BOOST_REQUIRE_EQUAL(svDatagroup.size(), 2u);
}

BOOST_AUTO_TEST_CASE( test_GetNodeRef)
{
    // test reference sequence of a segment

    const bam_header_info bamHeader(buildTestBamHeader());
    SVLocus locus;
    locus.addNode(GenomeInterval(0,1,1));
    GenomeInterval searchInterval;
    reference_contig_segment refSeq;
    getNodeRefSeq(bamHeader, locus, 0, getTestReferenceFilename(), searchInterval, refSeq);
    // check the size first
    BOOST_REQUIRE_EQUAL(refSeq.seq().size(), 101);
    // check the sequence
    BOOST_REQUIRE_EQUAL(refSeq.seq(), "GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGTATTTTCGTCTGGGGGGTGTGCACGCGATAGCATTGCGAGACGCTGG");
}

BOOST_AUTO_TEST_CASE( test_IsCandidateCountSufficient )
{
    // test candidates must have at least evidence of 2

    SVCandidate candidate;
    for (unsigned i(0); i<SVEvidenceType::SIZE; ++i)
        candidate.bp1.lowresEvidence.add(i,1);
    BOOST_REQUIRE(!isCandidateCountSufficient(candidate));
    for (unsigned i(0); i<SVEvidenceType::SIZE; ++i)
        candidate.bp1.lowresEvidence.add(i,1);
    BOOST_REQUIRE(isCandidateCountSufficient(candidate));


}

BOOST_AUTO_TEST_CASE( test_AddReadToDepthEst )
{
    // test depth

    bam_record bamRecord1;
    buildTestBamRecord(bamRecord1, 0, 200, 0, 210, 20, 15, "15M");
    bamRecord1.set_qname("Read-1");

    bam_record bamRecord2;
    buildTestBamRecord(bamRecord2, 0, 210, 0, 220, 20, 15, "15M");
    bamRecord2.set_qname("Read-2");

    std::vector<unsigned>depth(200);
    addReadToDepthEst(bamRecord1, 200, depth);
    addReadToDepthEst(bamRecord2, 200, depth);
    for (unsigned i = 0; i < 30; i++)
    {
        if (i >= 10 && i <= 19)
            BOOST_REQUIRE_EQUAL(depth[i], 2u); // second bamRead starts 10 bases after first bamRead
        else
            BOOST_REQUIRE_EQUAL(depth[i], 1u);
    }
}

BOOST_AUTO_TEST_CASE( test_IsBreakPointSignificant )
{
    // test the significance of a break point relative to a background noise rate

    std::vector<double> signalReadInfo;

    // minimum signal count should be 2
    BOOST_REQUIRE(!isBreakPointSignificant(0.1, 0.5, signalReadInfo));

    // Break point is not significant
    signalReadInfo.push_back(96);
    signalReadInfo.push_back(158);
    signalReadInfo.push_back(163);
    BOOST_REQUIRE(!isBreakPointSignificant(0.005, 0.00557491, signalReadInfo));

    // Breakpoint is significant
    signalReadInfo.clear();
    signalReadInfo.push_back(3440);
    signalReadInfo.push_back(3443);
    signalReadInfo.push_back(3452);
    signalReadInfo.push_back(3489);
    BOOST_REQUIRE(isBreakPointSignificant(0.03, 0.00836237, signalReadInfo));
}

BOOST_AUTO_TEST_CASE( test_IsSpanningCandidateSignalSignificant )
{
    // test the significance of a spanning candidate for  minimum supporting evidence

    SVCandidate svCandidate;
    FatSVCandidate fatSVCandidate(svCandidate, 1u);
    BOOST_REQUIRE(!isSpanningCandidateSignalSignificant(0.008, fatSVCandidate, 0));


    fatSVCandidate.bp1EvidenceIndex[0][0].push_back(3443);
    fatSVCandidate.bp1EvidenceIndex[0][0].push_back(3452);
    fatSVCandidate.bp1EvidenceIndex[0][0].push_back(3440);
    fatSVCandidate.bp1EvidenceIndex[0][0].push_back(3489);

    fatSVCandidate.bp2EvidenceIndex[0][0].push_back(1403);
    fatSVCandidate.bp2EvidenceIndex[0][0].push_back(1428);
    fatSVCandidate.bp2EvidenceIndex[0][0].push_back(1480);
    fatSVCandidate.bp2EvidenceIndex[0][0].push_back(1507);
    BOOST_REQUIRE(isSpanningCandidateSignalSignificant(0.008, fatSVCandidate, 0));
}

BOOST_AUTO_TEST_CASE( test_IsComplexCandidateSignalSignificant )
{
    // test the significance of a complex candidate for  minimum supporting evidence

    SVCandidate svCandidate;
    FatSVCandidate fatSVCandidate(svCandidate, 1u);
    fatSVCandidate.bp1EvidenceIndex[0][0].push_back(3443);
    fatSVCandidate.bp1EvidenceIndex[0][0].push_back(3452);
    fatSVCandidate.bp1EvidenceIndex[0][0].push_back(3440);
    fatSVCandidate.bp1EvidenceIndex[0][0].push_back(3489);
    BOOST_REQUIRE(!isComplexCandidateSignalSignificant(0.008, fatSVCandidate, 0));
}

BOOST_AUTO_TEST_CASE( test_IsAnySpanningCandidateSignalSignificant )
{
    // test the significance of spanning candidate across all the bams
    // relative to spanning noise rate.

    SVCandidate svCandidate;
    FatSVCandidate fatSVCandidate(svCandidate, 1u);

    fatSVCandidate.bp1EvidenceIndex[0][0].push_back(3443);
    fatSVCandidate.bp1EvidenceIndex[0][0].push_back(3452);
    fatSVCandidate.bp1EvidenceIndex[0][0].push_back(3440);
    fatSVCandidate.bp1EvidenceIndex[0][0].push_back(3489);

    fatSVCandidate.bp2EvidenceIndex[0][0].push_back(1403);
    fatSVCandidate.bp2EvidenceIndex[0][0].push_back(1428);
    fatSVCandidate.bp2EvidenceIndex[0][0].push_back(1480);
    fatSVCandidate.bp2EvidenceIndex[0][0].push_back(1507);
    std::vector<double > spanningNoiseRate;

    // spanning candidate is significant
    spanningNoiseRate.push_back(0.008);
    BOOST_REQUIRE(isAnySpanningCandidateSignalSignificant(1, fatSVCandidate, spanningNoiseRate));

    // not at all significant
    spanningNoiseRate.clear();
    spanningNoiseRate.push_back(0.1);
    BOOST_REQUIRE(!isAnySpanningCandidateSignalSignificant(1, fatSVCandidate, spanningNoiseRate));
}

BOOST_AUTO_TEST_CASE( test_IsAnyComplexCandidateSignalSignificant )
{
    // test the significance of complex candidate across all the bams
    // relative to assembly noise rate.

    SVCandidate svCandidate;
    FatSVCandidate fatSVCandidate(svCandidate, 1u);

    fatSVCandidate.bp1EvidenceIndex[0][0].push_back(3443);
    fatSVCandidate.bp1EvidenceIndex[0][0].push_back(3452);
    fatSVCandidate.bp1EvidenceIndex[0][0].push_back(3440);
    fatSVCandidate.bp1EvidenceIndex[0][0].push_back(3489);

    fatSVCandidate.bp2EvidenceIndex[0][0].push_back(1403);
    fatSVCandidate.bp2EvidenceIndex[0][0].push_back(1428);
    fatSVCandidate.bp2EvidenceIndex[0][0].push_back(1480);
    fatSVCandidate.bp2EvidenceIndex[0][0].push_back(1507);
    std::vector<double > assemblyNoiseRate;

    // complex candidate is significant
    assemblyNoiseRate.push_back(0.000002);
    BOOST_REQUIRE(isAnyComplexCandidateSignalSignificant(1, fatSVCandidate, assemblyNoiseRate));

    // complex candidate is not significant
    assemblyNoiseRate.clear();
    assemblyNoiseRate.push_back(0.008);
    BOOST_REQUIRE(!isAnyComplexCandidateSignalSignificant(1, fatSVCandidate, assemblyNoiseRate));
}

BOOST_AUTO_TEST_CASE( test_IsFilterSingleJunctionCandidate )
{
    // test the junction filter value

    SVCandidate svCandidate;
    FatSVCandidate fatSVCandidate(svCandidate, 1u);

    fatSVCandidate.bp1EvidenceIndex[0][0].push_back(3443);
    fatSVCandidate.bp1EvidenceIndex[0][0].push_back(3452);
    fatSVCandidate.bp1EvidenceIndex[0][0].push_back(3440);
    fatSVCandidate.bp1EvidenceIndex[0][0].push_back(3489);

    fatSVCandidate.bp2EvidenceIndex[0][0].push_back(1403);
    fatSVCandidate.bp2EvidenceIndex[0][0].push_back(1428);
    fatSVCandidate.bp2EvidenceIndex[0][0].push_back(1480);
    fatSVCandidate.bp2EvidenceIndex[0][0].push_back(1507);
    std::vector<double > assemblyNoiseRate;

    assemblyNoiseRate.push_back(0.000002);
    std::vector<double > spanningNoiseRate;
    spanningNoiseRate.push_back(0.008);

    // semi mapped
    BOOST_REQUIRE_EQUAL(isFilterSingleJunctionCandidate(false, spanningNoiseRate, assemblyNoiseRate, fatSVCandidate, 1), 1);

    svCandidate.bp1.state = SVBreakendState::index_t::RIGHT_OPEN ;
    svCandidate.bp2.state = SVBreakendState::index_t::LEFT_OPEN ;
    svCandidate.bp1.lowresEvidence.add(0, 2);
    FatSVCandidate fatSVCandidate1(svCandidate,1);

    // None
    BOOST_REQUIRE_EQUAL(isFilterSingleJunctionCandidate(true, spanningNoiseRate, assemblyNoiseRate, fatSVCandidate1, 1), 0);

    // spanning Low signal
    BOOST_REQUIRE_EQUAL(isFilterSingleJunctionCandidate(false, spanningNoiseRate, assemblyNoiseRate, fatSVCandidate1, 1), 4);

    // COMPLEX LOW COUNT
    svCandidate.bp1.state = SVBreakendState::index_t::COMPLEX ;
    svCandidate.bp2.state = SVBreakendState::index_t::UNKNOWN ;
    svCandidate.bp1.lowresEvidence.add(0, 2);
    FatSVCandidate fatSVCandidate2(svCandidate,1);
    BOOST_REQUIRE_EQUAL(isFilterSingleJunctionCandidate(false, spanningNoiseRate, assemblyNoiseRate, fatSVCandidate2, 1), 2);

    // COMPLEX LOW SIGNAL
    svCandidate.bp1.state = SVBreakendState::index_t::COMPLEX ;
    svCandidate.bp2.state = SVBreakendState::index_t::UNKNOWN ;
    svCandidate.bp1.lowresEvidence.clear();
    svCandidate.bp1.lowresEvidence.add(2, 3);
    FatSVCandidate fatSVCandidate3(svCandidate,1);
    BOOST_REQUIRE_EQUAL(isFilterSingleJunctionCandidate(false, spanningNoiseRate, assemblyNoiseRate, fatSVCandidate3, 1), 3);
}

BOOST_AUTO_TEST_CASE( test_filterCandidates )
{
    // test filters on all sv candidates

    SVCandidate svCandidate;
    FatSVCandidate fatSVCandidate(svCandidate, 1u);

    fatSVCandidate.bp1EvidenceIndex[0][0].push_back(3443);
    fatSVCandidate.bp1EvidenceIndex[0][0].push_back(3452);
    fatSVCandidate.bp1EvidenceIndex[0][0].push_back(3440);
    fatSVCandidate.bp1EvidenceIndex[0][0].push_back(3489);

    fatSVCandidate.bp2EvidenceIndex[0][0].push_back(1403);
    fatSVCandidate.bp2EvidenceIndex[0][0].push_back(1428);
    fatSVCandidate.bp2EvidenceIndex[0][0].push_back(1480);
    fatSVCandidate.bp2EvidenceIndex[0][0].push_back(1507);
    std::vector<double > assemblyNoiseRate;

    assemblyNoiseRate.push_back(0.000002);
    std::vector<double > spanningNoiseRate;
    spanningNoiseRate.push_back(0.008);

    std::vector<FatSVCandidate> svs;
    // semi mapped
    svs.push_back(fatSVCandidate);

    // spanning Low signal
    SVCandidate svCandidate1;
    svCandidate1.bp1.state = SVBreakendState::index_t::RIGHT_OPEN ;
    svCandidate1.bp2.state = SVBreakendState::index_t::LEFT_OPEN ;
    svCandidate1.bp1.lowresEvidence.add(0, 2);
    FatSVCandidate fatSVCandidate1(svCandidate1,1);
    svs.push_back(fatSVCandidate1);

    // COMPLEX LOW COUNT
    SVCandidate svCandidate2;
    svCandidate2.bp1.state = SVBreakendState::index_t::COMPLEX ;
    svCandidate2.bp2.state = SVBreakendState::index_t::UNKNOWN ;
    svCandidate2.bp1.lowresEvidence.add(0, 2);
    FatSVCandidate fatSVCandidate2(svCandidate2,1);
    svs.push_back(fatSVCandidate2);

    // COMPLEX LOW SIGNAL
    SVCandidate svCandidate3;
    svCandidate3.bp1.state = SVBreakendState::index_t::COMPLEX ;
    svCandidate3.bp2.state = SVBreakendState::index_t::UNKNOWN ;
    svCandidate3.bp1.lowresEvidence.clear();
    svCandidate3.bp1.lowresEvidence.add(2, 3);
    FatSVCandidate fatSVCandidate3(svCandidate3,1);
    svs.push_back(fatSVCandidate3);
    SVFinderStats stats;
    filterCandidates(false, spanningNoiseRate, assemblyNoiseRate, svs, stats, 1);

    // check all the stats
    BOOST_REQUIRE_EQUAL(svs.size(), 1);
    BOOST_REQUIRE_EQUAL(stats.ComplexLowCountFilter, 1);
    BOOST_REQUIRE_EQUAL(stats.ComplexLowSignalFilter, 1);
    BOOST_REQUIRE_EQUAL(stats.semiMappedFilter, 1);

    // check whether spanning Low signal sv candidate is there or not
    BOOST_REQUIRE_EQUAL(svs[0].bp1.state, SVBreakendState::index_t::RIGHT_OPEN);
    BOOST_REQUIRE_EQUAL(svs[0].bp2.state, SVBreakendState::index_t::LEFT_OPEN);
}

BOOST_AUTO_TEST_CASE( test_updateEvidenceIndex )
{
    // test the additional read info like read index

    bam_record bamRecord1;
    buildTestBamRecord(bamRecord1, 0, 200, 0, 210, 20, 15, "15M");
    bamRecord1.set_qname("Read-1");

    SVCandidateSetSequenceFragment fragment;
    SVObservation obs;
    obs.dnaFragmentSVEvidenceSource = SourceOfSVEvidenceInDNAFragment::READ1;
    fragment.read1.bamrec = bamRecord1;
    fragment.read1.readIndex = 1;
    obs.svEvidenceType = SVEvidenceType::SEMIALIGN;
    SVCandidate sv;
    FatSVCandidate fsv(sv, 1u);
    updateEvidenceIndex(fragment, obs, fsv, 0);
    BOOST_REQUIRE_EQUAL(fsv.bp1EvidenceIndex[SVEvidenceType::SEMIALIGN][0][0], 1);

    obs.svEvidenceType = SVEvidenceType::SPLIT_ALIGN;
    updateEvidenceIndex(fragment, obs, fsv, 0);
    BOOST_REQUIRE_EQUAL(fsv.bp1EvidenceIndex[SVEvidenceType::SPLIT_ALIGN][0].size(), 1);
    BOOST_REQUIRE_EQUAL(fsv.bp2EvidenceIndex[SVEvidenceType::SPLIT_ALIGN][0].size(), 0);
    BOOST_REQUIRE_EQUAL(fsv.bp1EvidenceIndex[SVEvidenceType::SPLIT_ALIGN][0][0], 1);

    // Adding supplementary read
    bam_record supplementSASplitRead;
    buildTestBamRecord(supplementSASplitRead);
    addSupplementaryAlignmentEvidence(supplementSASplitRead);
    SVCandidateSetRead read;
    read.bamrec = supplementSASplitRead;
    fragment.read1Supplemental.push_back(read);
    updateEvidenceIndex(fragment, obs, fsv, 0);
    BOOST_REQUIRE_EQUAL(fsv.bp1EvidenceIndex[SVEvidenceType::SPLIT_ALIGN][0].size(), 2);
    BOOST_REQUIRE_EQUAL(fsv.bp2EvidenceIndex[SVEvidenceType::SPLIT_ALIGN][0].size(), 1);
    BOOST_REQUIRE_EQUAL(fsv.bp1EvidenceIndex[SVEvidenceType::SPLIT_ALIGN][0][1], 1);
    BOOST_REQUIRE_EQUAL(fsv.bp2EvidenceIndex[SVEvidenceType::SPLIT_ALIGN][0][0], 0);

    // Read_PAIR
    bam_record bamRecord2;
    buildTestBamRecord(bamRecord2, 0, 210, 0, 220, 20, 15, "15M");
    bamRecord2.set_qname("Read-2");
    obs.dnaFragmentSVEvidenceSource = SourceOfSVEvidenceInDNAFragment::READ_PAIR;
    fragment.read1Supplemental.clear();
    fragment.read2.bamrec = bamRecord2;
    fragment.read2.readIndex = 2;
    obs.svEvidenceType = SVEvidenceType::PAIR;
    updateEvidenceIndex(fragment, obs, fsv, 0);
    BOOST_REQUIRE_EQUAL(fsv.bp1EvidenceIndex[SVEvidenceType::PAIR][0].size(), 1);
    BOOST_REQUIRE_EQUAL(fsv.bp2EvidenceIndex[SVEvidenceType::PAIR][0].size(), 1);
    BOOST_REQUIRE_EQUAL(fsv.bp1EvidenceIndex[SVEvidenceType::PAIR][0][0], 1);
    BOOST_REQUIRE_EQUAL(fsv.bp2EvidenceIndex[SVEvidenceType::PAIR][0][0], 2);
}

BOOST_AUTO_TEST_SUITE_END()