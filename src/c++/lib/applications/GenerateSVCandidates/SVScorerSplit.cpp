// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Manta
// Copyright (c) 2013 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/sequencing/licenses/>
//

///
/// \author Chris Saunders and Xiaoyu Chen
///

#include "SVScorer.hh"
#include "SVScorerShared.hh"


//#define DEBUG_SVS

#ifdef DEBUG_SVS
#include <iostream>
#include "blt_util/log.hh"
#endif



static
void
incrementAlleleEvidence(
    const SplitReadAlignment& bp1SR,
    const SplitReadAlignment& bp2SR,
    const unsigned readMapQ,
    SVSampleAlleleInfo& allele,
    SVFragmentEvidenceAlleleBreakendPerRead& bp1Support,
    SVFragmentEvidenceAlleleBreakendPerRead& bp2Support)
{
    float bp1Evidence(0);
    float bp2Evidence(0);
    if (bp1SR.has_evidence())
    {
        bp1Evidence = bp1SR.get_evidence();
        bp1Support.isSplitSupport = true;
        bp1Support.splitEvidence = bp1Evidence;
    }

    bp1Support.splitLnLhood = bp1SR.get_alignment().get_alignLnLhood();

    if (bp2SR.has_evidence())
    {
        bp2Evidence = bp2SR.get_evidence();
        bp2Support.isSplitSupport = true;
        bp2Support.splitEvidence = bp2Evidence;
    }

    bp2Support.splitLnLhood = bp2SR.get_alignment().get_alignLnLhood();

    const float evidence(std::max(bp1Evidence, bp2Evidence));

    if ((bp1SR.has_evidence()) || (bp2SR.has_evidence()))
    {
        allele.splitReadCount++;
        allele.splitReadEvidence += evidence;
        allele.splitReadMapQ += readMapQ * readMapQ;

#ifdef DEBUG_SVS
        log_os << "bp1\n";
        log_os << bp1SR;
        log_os << "bp2\n";
        log_os << bp2SR;
        log_os << "evidence = " << evidence << "\n";
        log_os << "accumulated evidence = " << allele.splitReadEvidence << "\n";
        log_os << "contigCount = " << allele.splitReadCount << "\n\n";
#endif
    }
}



static
void
scoreSplitReads(
    const CallOptionsSharedDeriv& dopt,
    const SVBreakend& bp,
    const SVAlignmentInfo& svAlignInfo,
    const unsigned minMapQ,
    SVEvidence::evidenceTrack_t& sampleEvidence,
    bam_streamer& readStream,
    SVSampleInfo& sample)
{
    // extract reads overlapping the break point
    readStream.set_new_region(bp.interval.tid, bp.interval.range.begin_pos(), bp.interval.range.end_pos());
    while (readStream.next())
    {
        const bam_record& bamRead(*(readStream.get_record_ptr()));

        if (bamRead.is_filter()) continue;
        if (bamRead.is_dup()) continue;
        if (bamRead.is_secondary()) continue;
        if (bamRead.is_supplement()) continue;

        const std::string readSeq = bamRead.get_bam_read().get_string();
        const uint8_t* qual(bamRead.qual());
        const unsigned readMapQ = bamRead.map_qual();

        SVFragmentEvidence& fragment(sampleEvidence[bamRead.qname()]);

        const bool isRead1(bamRead.is_first());
        setReadEvidence(minMapQ, bamRead, fragment.getRead(isRead1));

        SVFragmentEvidenceAlleleBreakendPerRead& altBp1ReadSupport(fragment.alt.bp1.getRead(isRead1));
        SVFragmentEvidenceAlleleBreakendPerRead& refBp1ReadSupport(fragment.ref.bp1.getRead(isRead1));
        SVFragmentEvidenceAlleleBreakendPerRead& altBp2ReadSupport(fragment.alt.bp2.getRead(isRead1));
        SVFragmentEvidenceAlleleBreakendPerRead& refBp2ReadSupport(fragment.ref.bp2.getRead(isRead1));

        /// in this function we evaluate the hypothesis of both breakends at the same time, the only difference bp1 vs
        /// bp2 makes is where in the bam we look for reads, therefore if we see split evaluation for bp1 or bp2, we can skip this read:
        if (altBp1ReadSupport.isSplitEvaluated) continue;

        altBp1ReadSupport.isSplitEvaluated = true;
        refBp1ReadSupport.isSplitEvaluated = true;
        altBp2ReadSupport.isSplitEvaluated = true;
        refBp2ReadSupport.isSplitEvaluated = true;

        // align the read to the somatic contig
        SplitReadAlignment bp1ContigSR;
        bp1ContigSR.align(readSeq, dopt.altQ, qual, svAlignInfo.bp1ContigSeq(), svAlignInfo.bp1ContigOffset);
        SplitReadAlignment bp2ContigSR;
        bp2ContigSR.align(readSeq, dopt.altQ, qual, svAlignInfo.bp2ContigSeq(), svAlignInfo.bp2ContigOffset);

        // align the read to reference regions
        SplitReadAlignment bp1RefSR;
        bp1RefSR.align(readSeq, dopt.refQ, qual, svAlignInfo.bp1ReferenceSeq(), svAlignInfo.bp1RefOffset);
        SplitReadAlignment bp2RefSR;
        bp2RefSR.align(readSeq, dopt.refQ, qual, svAlignInfo.bp2ReferenceSeq(), svAlignInfo.bp2RefOffset);

        // scoring
        incrementAlleleEvidence(bp1ContigSR, bp2ContigSR, readMapQ, sample.alt, altBp1ReadSupport, altBp2ReadSupport);
        incrementAlleleEvidence(bp1RefSR, bp2RefSR, readMapQ, sample.ref, refBp1ReadSupport, refBp2ReadSupport);
    }
}



/// return rms given sum of squares
static
float
finishRms(
    const float sumSqr,
    const unsigned count)
{
    if (count == 0) return 0.;
    return std::sqrt(sumSqr / static_cast<float>(count));
}



static
void
finishRms(
    SVSampleAlleleInfo& sai)
{
    sai.splitReadMapQ = finishRms(sai.splitReadMapQ, sai.splitReadCount);
}



/// make final split read computations after bam scanning is finished:
static
void
finishSampleSRData(
    SVSampleInfo& sample)
{
    // finish rms mapq:
    finishRms(sample.alt);
    finishRms(sample.ref);
}



void
SVScorer::
getSVSplitReadSupport(
    const SVCandidateAssemblyData& assemblyData,
    const SVCandidate& sv,
    SVScoreInfo& baseInfo,
    SVEvidence& evidence)
{
    static const unsigned maxDepthSRFactor(2); ///< at what multiple of the maxDepth do we skip split read analysis?

    bool isSkipSRSearchDepth(false);

    if (_dFilterDiploid.isMaxDepthFilter() && _dFilterSomatic.isMaxDepthFilter())
    {
        const double bp1MaxMaxDepth(std::max(_dFilterDiploid.maxDepth(sv.bp1.interval.tid), _dFilterSomatic.maxDepth(sv.bp1.interval.tid)));
        const double bp2MaxMaxDepth(std::max(_dFilterDiploid.maxDepth(sv.bp2.interval.tid), _dFilterSomatic.maxDepth(sv.bp2.interval.tid)));

        isSkipSRSearchDepth=((baseInfo.bp1MaxDepth > (maxDepthSRFactor*bp1MaxMaxDepth)) ||
                             (baseInfo.bp2MaxDepth > (maxDepthSRFactor*bp2MaxMaxDepth)));
    }

    // apply the split-read scoring, only when:
    // 1) the SV is precise, i.e. has successful somatic contigs;
    // 2) the values of max depth are reasonable (otherwise, the read map may blow out).
    const bool isSkipSRSearch(
        (sv.isImprecise()) ||
        (isSkipSRSearchDepth));

    if (isSkipSRSearch) return;

    // Get Data on standard read pairs crossing the two breakends,

    // extract SV alignment info for split read evidence
    const SVAlignmentInfo SVAlignInfo(sv, assemblyData);
#ifdef DEBUG_SVS
    log_os << SVAlignInfo << "\n";
#endif

    const unsigned minMapQ(_readScanner.getMinMapQ());

    const unsigned bamCount(_bamStreams.size());
    for (unsigned bamIndex(0); bamIndex < bamCount; ++bamIndex)
    {
        const bool isTumor(_isAlignmentTumor[bamIndex]);
        SVSampleInfo& sample(isTumor ? baseInfo.tumor : baseInfo.normal);
        bam_streamer& bamStream(*_bamStreams[bamIndex]);

        SVEvidence::evidenceTrack_t& sampleEvidence(evidence.getSample(isTumor));

        // scoring split reads overlapping bp1
        scoreSplitReads(_callDopt, sv.bp1, SVAlignInfo, minMapQ, sampleEvidence,
                        bamStream, sample);
        // scoring split reads overlapping bp2
        scoreSplitReads(_callDopt, sv.bp2, SVAlignInfo, minMapQ, sampleEvidence,
                        bamStream, sample);
    }

    finishSampleSRData(baseInfo.tumor);
    finishSampleSRData(baseInfo.normal);

#ifdef DEBUG_SVS
    log_os << "tumor contig SP count: " << baseInfo.tumor.alt.splitReadCount << "\n";
    log_os << "tumor contig SP evidence: " << baseInfo.tumor.alt.splitReadEvidence << "\n";
    log_os << "tumor contig SP_mapQ: " << baseInfo.tumor.alt.splitReadMapQ << "\n";
    log_os << "normal contig SP count: " << baseInfo.normal.alt.splitReadCount << "\n";
    log_os << "normal contig SP evidence: " << baseInfo.normal.alt.splitReadEvidence << "\n";
    log_os << "normal contig SP_mapQ: " << baseInfo.normal.alt.splitReadMapQ << "\n";

    log_os << "tumor ref SP count: " << baseInfo.tumor.ref.splitReadCount << "\n";
    log_os << "tumor ref SP evidence: " << baseInfo.tumor.ref.splitReadEvidence << "\n";
    log_os << "tumor ref SP_mapQ: " << baseInfo.tumor.ref.splitReadMapQ << "\n";
    log_os << "normal ref SP count: " << baseInfo.normal.ref.splitReadCount << "\n";
    log_os << "normal ref SP evidence: " << baseInfo.normal.ref.splitReadEvidence << "\n";
    log_os << "normal ref SP_mapQ: " << baseInfo.normal.ref.splitReadMapQ << "\n";
#endif
}
