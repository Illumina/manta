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

#include "blt_util/align_path_bam_util.hh"
#include "blt_util/bam_streamer.hh"
#include "blt_util/prob_util.hh"
#include "blt_util/qscore.hh"
#include "common/Exceptions.hh"
#include "manta/ReadGroupStatsSet.hh"

#include "boost/foreach.hpp"

#include <algorithm>
#include <iostream>
#include <string>

//#define DEBUG_SVS

#ifdef DEBUG_SVS
#include "blt_util/log.hh"
#endif


//#define DEBUG_SCORE

#ifdef DEBUG_SCORE
#include "blt_util/log.hh"
#endif


const pos_t PairOptions::minFragSupport(50);



SVScorer::
SVScorer(
    const GSCOptions& opt,
    const bam_header_info& header) :
    _isAlignmentTumor(opt.alignFileOpt.isAlignmentTumor),
    _diploidOpt(opt.diploidOpt),
    _diploidDopt(_diploidOpt),
    _somaticOpt(opt.somaticOpt),
    _dFilterDiploid(opt.chromDepthFilename, _diploidOpt.maxDepthFactor, header),
    _dFilterSomatic(opt.chromDepthFilename, _somaticOpt.maxDepthFactor, header),
    _readScanner(opt.scanOpt,opt.statsFilename,opt.alignFileOpt.alignmentFilename)
{
    // setup regionless bam_streams:
    // setup all data for main analysis loop:
    BOOST_FOREACH(const std::string& afile, opt.alignFileOpt.alignmentFilename)
    {
        // avoid creating shared_ptr temporaries:
        streamPtr tmp(new bam_streamer(afile.c_str()));
        _bamStreams.push_back(tmp);
    }
}



/// add bam alignment to simple short-range vector depth estimate
///
/// \param[in] beginPos this is the begin position of the range covered by the depth array
///
static
void
addReadToDepthEst(
    const bam_record& bamRead,
    const pos_t beginPos,
    std::vector<unsigned>& depth)
{
    using namespace ALIGNPATH;

    const pos_t endPos(beginPos+depth.size());

    // get cigar:
    path_t apath;
    bam_cigar_to_apath(bamRead.raw_cigar(), bamRead.n_cigar(), apath);

    pos_t refPos(bamRead.pos()-1);
    BOOST_FOREACH(const path_segment& ps, apath)
    {
        if (refPos>=endPos) return;

        if (MATCH == ps.type)
        {
            for (pos_t pos(refPos); pos < (refPos+static_cast<pos_t>(ps.length)); ++pos)
            {
                if (pos>=beginPos)
                {
                    if (pos>=endPos) return;
                    depth[pos-beginPos]++;
                }
            }
        }
        if (is_segment_type_ref_length(ps.type)) refPos += ps.length;
    }
}



unsigned
SVScorer::
getBreakendMaxMappedDepth(
    const SVBreakend& bp)
{
    /// define a new interval -/+ 50 bases around the center pos
    /// of the breakpoint
    static const pos_t regionSize(50);
    const pos_t centerPos(bp.interval.range.center_pos());
    const known_pos_range2 searchRange(std::max((centerPos-regionSize),0), (centerPos+regionSize));

    std::vector<unsigned> depth(searchRange.size(),0);

    bool isNormalFound(false);

    const unsigned bamCount(_bamStreams.size());
    for (unsigned bamIndex(0); bamIndex < bamCount; ++bamIndex)
    {
        if (_isAlignmentTumor[bamIndex]) continue;

        bam_streamer& bamStream(*_bamStreams[bamIndex]);

        // set bam stream to new search interval:
        bamStream.set_new_region(bp.interval.tid, searchRange.begin_pos(), searchRange.end_pos());

        while (bamStream.next())
        {
            const bam_record& bamRead(*(bamStream.get_record_ptr()));

            // turn filtration off down to mapped only to match depth estimate method:
            //if (_readScanner.isReadFiltered(bamRead)) continue;
            if (bamRead.is_unmapped()) continue;

            if ((bamRead.pos()-1) >= searchRange.end_pos()) break;

            addReadToDepthEst(bamRead,searchRange.begin_pos(),depth);
        }

        isNormalFound=true;
        break;
    }

    assert(isNormalFound);

    return *(std::max_element(depth.begin(),depth.end()));
}



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

        setReadEvidence(minMapQ, bamRead, fragment.getRead(bamRead.is_first()));

        SVFragmentEvidenceAlleleBreakendPerRead& altBp1ReadSupport(fragment.alt.bp1.getRead(bamRead.is_first()));
        SVFragmentEvidenceAlleleBreakendPerRead& refBp1ReadSupport(fragment.ref.bp1.getRead(bamRead.is_first()));
        SVFragmentEvidenceAlleleBreakendPerRead& altBp2ReadSupport(fragment.alt.bp2.getRead(bamRead.is_first()));
        SVFragmentEvidenceAlleleBreakendPerRead& refBp2ReadSupport(fragment.ref.bp2.getRead(bamRead.is_first()));

        /// in this function we evaluate the hypothesis of both breakends at the same time, the only difference bp1 vs
        /// bp2 makes is where in the bam we look for reads, therefore if we see split evaluation for bp1 or bp2, we can skip this read:
        if (altBp1ReadSupport.isSplitEvaluated) continue;

        altBp1ReadSupport.isSplitEvaluated = true;
        refBp1ReadSupport.isSplitEvaluated = true;
        altBp2ReadSupport.isSplitEvaluated = true;
        refBp2ReadSupport.isSplitEvaluated = true;

        // align the read to the somatic contig
        SplitReadAlignment bp1ContigSR;
        bp1ContigSR.align(readSeq, qual, svAlignInfo.bp1ContigSeq(), svAlignInfo.bp1ContigOffset);
        SplitReadAlignment bp2ContigSR;
        bp2ContigSR.align(readSeq, qual, svAlignInfo.bp2ContigSeq(), svAlignInfo.bp2ContigOffset);

        // align the read to reference regions
        SplitReadAlignment bp1RefSR;
        bp1RefSR.align(readSeq, qual, svAlignInfo.bp1ReferenceSeq(), svAlignInfo.bp1RefOffset);
        SplitReadAlignment bp2RefSR;
        bp2RefSR.align(readSeq, qual, svAlignInfo.bp2ReferenceSeq(), svAlignInfo.bp2RefOffset);

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
        scoreSplitReads(sv.bp1, SVAlignInfo, minMapQ, sampleEvidence,
                        bamStream, sample);
        // scoring split reads overlapping bp2
        scoreSplitReads(sv.bp2, SVAlignInfo, minMapQ, sampleEvidence,
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



static
bool
isAlleleReadSupport(
    const SVFragmentEvidenceAllele& allele,
    const bool isRead1)
{
    return (allele.bp1.getRead(isRead1).isSplitSupport ||
            allele.bp2.getRead(isRead1).isSplitSupport);
}



static
bool
isAnyReadSupport(
    const SVFragmentEvidence& fragev,
    const bool isRead1)
{
    return (isAlleleReadSupport(fragev.alt, isRead1) ||
            isAlleleReadSupport(fragev.ref, isRead1));
}



static
void
lnToProb(
    float& lower,
    float& higher)
{
    lower = std::exp(lower-higher);
    higher = 1./(lower+1.);
    lower  = lower/(lower+1.);
}



static
void
addReadSupport(
    const SVFragmentEvidence& fragev,
    const bool isRead1,
    SVSampleInfo& sampleBaseInfo)
{
    static const float supportProb(0.9999);

    // only consider reads where at least one allele and one breakend is confident
    //
    // ...note this is done in the absence of having a noise state in the model
    //
    if (! isAnyReadSupport(fragev,isRead1)) return;

    float altLhood =
        std::max(fragev.alt.bp1.getRead(isRead1).splitLnLhood,
                 fragev.alt.bp2.getRead(isRead1).splitLnLhood);

    float refLhood =
        std::max(fragev.ref.bp1.getRead(isRead1).splitLnLhood,
                 fragev.ref.bp2.getRead(isRead1).splitLnLhood);

    // convert to normalized prob:
    if (altLhood > refLhood)
    {
        lnToProb(refLhood, altLhood);
        if (altLhood>supportProb) sampleBaseInfo.alt.confidentSplitReadCount++;
    }
    else
    {
        lnToProb(altLhood, refLhood);
        if (refLhood>supportProb) sampleBaseInfo.ref.confidentSplitReadCount++;
    }


}



static
void
getSampleCounts(
    const SVEvidence::evidenceTrack_t& sampleEvidence,
    SVSampleInfo& sampleBaseInfo)
{
    BOOST_FOREACH(const SVEvidence::evidenceTrack_t::value_type& val, sampleEvidence)
    {

        const SVFragmentEvidence& fragev(val.second);

        // evaluate read1 and read2 from this fragment
        //
        addReadSupport(fragev,true,sampleBaseInfo);
        addReadSupport(fragev,false,sampleBaseInfo);
    }
}



static
void
getSVSupportSummary(
    const SVEvidence& evidence,
    SVScoreInfo& baseInfo)
{
    /// get conservative count of reads which support only one allele, ie. P ( allele | read ) is high
    ///
    getSampleCounts(evidence.normal, baseInfo.normal);
    getSampleCounts(evidence.tumor, baseInfo.tumor);
}



/// shared information gathering steps of all scoring models
void
SVScorer::
scoreSV(
    const SVCandidateSetData& svData,
    const SVCandidateAssemblyData& assemblyData,
    const SVCandidate& sv,
    SVScoreInfo& baseInfo,
    SVEvidence& evidence)
{
    // get breakend center_pos depth estimate:
    baseInfo.bp1MaxDepth=(getBreakendMaxMappedDepth(sv.bp1));
    baseInfo.bp2MaxDepth=(getBreakendMaxMappedDepth(sv.bp2));

    /// global evidence accumulator for this SV:

    // count the paired-read fragments supporting the ref and alt alleles in each sample:
    //
    getSVPairSupport(svData, sv, baseInfo, evidence);

    // count the split reads supporting the ref and alt alleles in each sample
    //
    getSVSplitReadSupport(assemblyData, sv, baseInfo, evidence);

    // compute allele likelihoods, and any other summary metric shared between all models:
    //
    getSVSupportSummary(evidence, baseInfo);
}



static
void
incrementAlleleLhood(
    const float chimeraRate,
    const float chimeraComp,
    const SVFragmentEvidenceAllele& allele,
    float& bpLhood)
{
    float fragProb(0);
    if (allele.bp1.isFragmentSupport)
    {
        fragProb = allele.bp1.fragLengthProb;
    }

    if (allele.bp2.isFragmentSupport)
    {
        fragProb = std::max(fragProb, allele.bp2.fragLengthProb);
    }

    bpLhood *= (chimeraComp*fragProb + chimeraRate);
}



static
bool
isAnyPairSupport(
    const SVFragmentEvidence& fragev)
{
    const bool isRefSupport(fragev.ref.bp1.isFragmentSupport || fragev.ref.bp2.isFragmentSupport);
    const bool isAltSupport(fragev.alt.bp1.isFragmentSupport || fragev.alt.bp2.isFragmentSupport);

    return (isRefSupport || isAltSupport);
}



/// score diploid germline specific components:
static
void
scoreDiploidSV(
    const CallOptionsDiploid& diploidOpt,
    const CallOptionsDiploidDeriv& diploidDopt,
    const SVCandidate& sv,
    const ChromDepthFilterUtil& dFilter,
    const SVEvidence& evidence,
    SVScoreInfo& baseInfo,
    SVScoreInfoDiploid& diploidInfo)
{
#ifdef DEBUG_SCORE
    static const std::string logtag("scoreDiploidSV: ");
#endif

    /// TODO: set this from graph data:
    ///
    /// put some more thought into this -- is this P (spurious | any old read) or P( spurious | chimera ) ??
    /// it seems like it should be the later in the usages that really matter.
    ///
    static const float chimeraRate(1e-3);
    static const float chimeraComp(1.-chimeraRate);

    //
    // compute qualities
    //
    {
        float loglhood[DIPLOID_GT::SIZE];
        for (unsigned gt(0); gt<DIPLOID_GT::SIZE; ++gt)
        {
            loglhood[gt] = 0.;
        }

        BOOST_FOREACH(const SVEvidence::evidenceTrack_t::value_type& val, evidence.normal)
        {
            const SVFragmentEvidence& fragev(val.second);

            float fragRefLhood(1);
            float fragAltLhood(1);

#ifdef DEBUG_SCORE
            log_os << logtag << "qname: " << val.first << " fragev: " << fragev << "\n";
#endif

            /// high-quality spanning support relies on read1 and read2 mapping well:
            if ( fragev.read1.isObservedAnchor() && fragev.read2.isObservedAnchor())
            {
                /// only add to the likelihood if the fragment "supports" at least one allele:
                if ( isAnyPairSupport(fragev) )
                {
                    incrementAlleleLhood(chimeraRate, chimeraComp, fragev.alt, fragAltLhood);
                    incrementAlleleLhood(chimeraRate, chimeraComp, fragev.ref, fragRefLhood);
                }
            }

//            float read1RefLhood(1);
//            float read1AltLhood(1);

            for (unsigned gt(0); gt<DIPLOID_GT::SIZE; ++gt)
            {
                const float altFrac(DIPLOID_GT::altFraction(gt));
                const float reflhood(fragRefLhood * (1.-altFrac));
                const float altlhood(fragAltLhood * altFrac);
                loglhood[gt] += std::log(reflhood + altlhood);

#ifdef DEBUG_SCORE
                log_os << logtag << "gt/fragref/ref/fragalt/alt: "
                       << DIPLOID_GT::label(gt)
                       << " " << fragRefLhood
                       << " " << reflhood
                       << " " << fragAltLhood
                       << " " << altlhood
                       << "\n";
#endif
            }
        }

        float pprob[DIPLOID_GT::SIZE];
        for (unsigned gt(0); gt<DIPLOID_GT::SIZE; ++gt)
        {
            pprob[gt] = loglhood[gt] + diploidDopt.prior[gt];
        }

        unsigned maxGt(0);
        normalize_ln_distro(pprob, pprob+DIPLOID_GT::SIZE, maxGt);

#ifdef DEBUG_SCORE
        for (unsigned gt(0); gt<DIPLOID_GT::SIZE; ++gt)
        {
            log_os << logtag << "gt/lhood/prior/pprob: "
                   << DIPLOID_GT::label(gt)
                   << " " << loglhood[gt]
                   << " " << diploidDopt.prior[gt]
                   << " " << pprob[gt]
                   << "\n";
        }
#endif

        diploidInfo.gt=static_cast<DIPLOID_GT::index_t>(maxGt);
        diploidInfo.altScore=error_prob_to_qphred(pprob[DIPLOID_GT::REF]);
        diploidInfo.gtScore=error_prob_to_qphred(prob_comp(pprob,pprob+DIPLOID_GT::SIZE, diploidInfo.gt));
    }


    //
    // apply filters
    //
    if (diploidInfo.altScore >= diploidOpt.minOutputAltScore)
    {
        if (dFilter.isMaxDepthFilter())
        {
            // apply maxdepth filter if either of the breakpoints exceeds the maximum depth:
            if (baseInfo.bp1MaxDepth > dFilter.maxDepth(sv.bp1.interval.tid))
            {
                diploidInfo.filters.insert(diploidOpt.maxDepthFilterLabel);
            }
            else if (baseInfo.bp2MaxDepth > dFilter.maxDepth(sv.bp2.interval.tid))
            {
                diploidInfo.filters.insert(diploidOpt.maxDepthFilterLabel);
            }
        }

        if ( diploidInfo.gtScore < diploidOpt.minGTScoreFilter)
        {
            diploidInfo.filters.insert(diploidOpt.minGTFilterLabel);
        }
    }
}



/// score somatic specific components:
static
void
scoreSomaticSV(
    const CallOptionsSomatic& somaticOpt,
    const SVCandidate& sv,
    const ChromDepthFilterUtil& dFilter,
    SVScoreInfo& baseInfo,
    SVScoreInfoSomatic& somaticInfo)
{
    //
    // compute qualities
    //
    {
        bool isNonzeroSomaticQuality(true);

        /// first check for substantial support in the normal:
        if (baseInfo.normal.alt.spanPairCount > 1) isNonzeroSomaticQuality=false;
        if (baseInfo.normal.alt.confidentSplitReadCount > 5) isNonzeroSomaticQuality=false;

        if (isNonzeroSomaticQuality)
        {
            const bool lowPairSupport(baseInfo.tumor.alt.spanPairCount < 6);
            const bool lowSplitSupport(baseInfo.tumor.alt.confidentSplitReadCount < 6);
            const bool lowSingleSupport((baseInfo.tumor.alt.bp1SpanReadCount < 14) || (baseInfo.tumor.alt.bp2SpanReadCount < 14));
            const bool highSingleContam((baseInfo.normal.alt.bp1SpanReadCount > 1) || (baseInfo.normal.alt.bp2SpanReadCount > 1));

            /// allow single pair support to rescue an SV only if the evidence looks REALLY good:
            if ((lowPairSupport && lowSplitSupport) && (lowSingleSupport || highSingleContam))
                isNonzeroSomaticQuality=false;
        }

        if (isNonzeroSomaticQuality)
        {
            if (baseInfo.normal.alt.spanPairCount)
            {
                const double ratio(static_cast<double>(baseInfo.tumor.alt.spanPairCount)/static_cast<double>(baseInfo.normal.alt.spanPairCount));
                if (ratio<9)
                {
                    isNonzeroSomaticQuality=false;
                }
            }
            if (baseInfo.normal.alt.bp1SpanReadCount)
            {
                const double ratio(static_cast<double>(baseInfo.tumor.alt.bp1SpanReadCount)/static_cast<double>(baseInfo.normal.alt.bp1SpanReadCount));
                if (ratio<9)
                {
                    isNonzeroSomaticQuality=false;
                }
            }
            if (baseInfo.normal.alt.bp2SpanReadCount)
            {
                const double ratio(static_cast<double>(baseInfo.tumor.alt.bp2SpanReadCount)/static_cast<double>(baseInfo.normal.alt.bp2SpanReadCount));
                if (ratio<9)
                {
                    isNonzeroSomaticQuality=false;
                }
            }
        }

        {
            // there needs to be some ref support in the normal as well:
            const bool normRefPairSupport(baseInfo.normal.ref.spanPairCount > 6);
            const bool normRefSplitSupport(baseInfo.normal.ref.confidentSplitReadCount > 6);

            if (! (normRefPairSupport || normRefSplitSupport)) isNonzeroSomaticQuality=false;
        }

        if (isNonzeroSomaticQuality) somaticInfo.somaticScore=60;
    }


    //
    // apply filters
    //
    {
        if (dFilter.isMaxDepthFilter())
        {
            // apply maxdepth filter if either of the breakpoints exceeds the maximum depth:
            if (baseInfo.bp1MaxDepth > dFilter.maxDepth(sv.bp1.interval.tid))
            {
                somaticInfo.filters.insert(somaticOpt.maxDepthFilterLabel);
            }
            else if (baseInfo.bp2MaxDepth > dFilter.maxDepth(sv.bp2.interval.tid))
            {
                somaticInfo.filters.insert(somaticOpt.maxDepthFilterLabel);
            }
        }
    }
}



void
SVScorer::
scoreSV(
    const SVCandidateSetData& svData,
    const SVCandidateAssemblyData& assemblyData,
    const SVCandidate& sv,
    const bool isSomatic,
    SVModelScoreInfo& modelScoreInfo)
{
    modelScoreInfo.clear();

    // accumulate model-neutral evidence for each candidate (or its corresponding reference allele)
    SVEvidence evidence;
    scoreSV(svData, assemblyData, sv, modelScoreInfo.base, evidence);

    // score components specific to diploid-germline model:
    scoreDiploidSV(_diploidOpt, _diploidDopt, sv, _dFilterDiploid, evidence, modelScoreInfo.base, modelScoreInfo.diploid);

    // score components specific to somatic model:
    if (isSomatic)
    {
        scoreSomaticSV(_somaticOpt, sv, _dFilterSomatic, modelScoreInfo.base, modelScoreInfo.somatic);
    }
}




