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

#include "blt_util/align_path_bam_util.hh"
#include "blt_util/bam_streamer.hh"
#include "blt_util/bam_record_util.hh"
#include "blt_util/log.hh"
#include "blt_util/ReadKey.hh"
#include "common/Exceptions.hh"
#include "manta/ReadGroupStatsSet.hh"
#include "boost/foreach.hpp"

#include <algorithm>
#include <iostream>
#include <string>

//#define DEBUG_SVS


SVScorer::
SVScorer(
    const GSCOptions& opt,
    const bam_header_info& header) :
    _isAlignmentTumor(opt.alignFileOpt.isAlignmentTumor),
    _diploidOpt(opt.diploidOpt),
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
    if (bp2SR.has_evidence())
    {
        bp2Evidence = bp2SR.get_evidence();
        bp2Support.isSplitSupport = true;
        bp2Support.splitEvidence = bp2Evidence;
    }

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
setReadEvidence(
    const unsigned minMapQ,
    const bam_record& bamRead,
    SVFragmentEvidenceRead& read)
{
    if (read.isScanned) return;

    read.isScanned = true;
    read.mapq = bamRead.map_qual();
    read.isAnchored = (read.mapq >= minMapQ);
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

        const std::string readSeq = bamRead.get_bam_read().get_string();
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
        bp1ContigSR.align(readSeq, svAlignInfo.bp1ContigSeq(), svAlignInfo.bp1ContigOffset);
        SplitReadAlignment bp2ContigSR;
        bp2ContigSR.align(readSeq, svAlignInfo.bp2ContigSeq(), svAlignInfo.bp2ContigOffset);

        // align the read to reference regions
        SplitReadAlignment bp1RefSR;
        bp1RefSR.align(readSeq, svAlignInfo.bp1RefSeq, svAlignInfo.bp1RefOffset);
        SplitReadAlignment bp2RefSR;
        bp2RefSR.align(readSeq, svAlignInfo.bp2RefSeq, svAlignInfo.bp2RefOffset);

        // scoring
        incrementAlleleEvidence(bp1ContigSR, bp2ContigSR, readMapQ, sample.alt, altBp1ReadSupport, altBp2ReadSupport);
        incrementAlleleEvidence(bp1RefSR, bp2RefSR, readMapQ, sample.ref, refBp1ReadSupport, refBp2ReadSupport);
    }
}



/// get reference allele support at a single breakend:
///
void
SVScorer::
getSVRefPairSupport(
    const SVBreakend& bp,
    SVScoreInfo& baseInfo,
    SVEvidence& evidence,
    const bool isBp1)
{
    /// search for all read pairs supporting the reference allele
    ///
    /// APPROXIMATION: for imprecise and precise variants treat the breakend locations as the center of the
    ///  breakend interval.
    ///
    /// TODO: improve on the approx above
    ///

    /// TODO: track read key to account for overlap of spanning and split read evidence

    const pos_t centerPos(bp.interval.range.center_pos());


    /// we're interested in any fragments which cross center pos with at least N bases of support on each side
    /// (note this definition is certain to overlap the split read definition whenever N is less than the read length

    static const pos_t minFragSupport(80);
    const unsigned minMapQ(_readScanner.getMinMapQ());

    const unsigned bamCount(_bamStreams.size());
    for (unsigned bamIndex(0); bamIndex < bamCount; ++bamIndex)
    {
        const bool isTumor(_isAlignmentTumor[bamIndex]);
        SVSampleInfo& sample(isTumor ? baseInfo.tumor : baseInfo.normal);

        bam_streamer& bamStream(*_bamStreams[bamIndex]);

        /// set the search range around centerPos so that we can get any fragments at the Xth percentile length or smaller which could have
        /// min Fragsupport
        const SVLocusScanner::Range& pRange(_readScanner.getEvidencePairRange(bamIndex));
        const unsigned minFrag(pRange.min);
        const unsigned maxFrag(pRange.max);

        const uint32_t beginPos(centerPos-maxFrag+minFragSupport);
        const uint32_t endPos(centerPos+maxFrag-minFragSupport+1);

        if (beginPos >= endPos) continue;

        // set bam stream to new search interval:
        bamStream.set_new_region(bp.interval.tid, beginPos, endPos);

        while (bamStream.next())
        {
            const bam_record& bamRead(*(bamStream.get_record_ptr()));

            if (bamRead.is_filter()) continue;
            if (bamRead.is_dup()) continue;
            if (bamRead.is_secondary()) continue;

            if (bamRead.is_unmapped() || bamRead.is_mate_unmapped()) continue;

            /// check for standard innie orientation:
            if (! is_innie_pair(bamRead)) continue;

            /// check if fragment is too big or too small:
            const unsigned tSize(std::abs(bamRead.template_size()));
            if (tSize < minFrag) continue;
            if (tSize > maxFrag) continue;

            // count only from the down stream read unless the mate-pos goes past center-pos
            const bool isLeftMost(bamRead.pos() < bamRead.mate_pos());
            const bool isRead1Tie((bamRead.pos() == bamRead.mate_pos()) && bamRead.is_first());
            const bool isDefaultSelected(isLeftMost || isRead1Tie);

            const bool isMateBeforeCenter(bamRead.mate_pos() < centerPos);

            bool isDoubleCountSkip(false);
            if ( isDefaultSelected && isMateBeforeCenter ) isDoubleCountSkip=true;
            if ( (!isDefaultSelected) && (!isMateBeforeCenter) ) isDoubleCountSkip=true;

            // get fragment range:
            pos_t fragBegin(0);
            if (isLeftMost)
            {
                fragBegin=bamRead.pos()-1;
            }
            else
            {
                fragBegin=bamRead.mate_pos()-1;
            }
            const pos_t fragEnd(fragBegin+std::abs(bamRead.template_size()));

            if (fragBegin > fragEnd)
            {
                using namespace illumina::common;

                std::ostringstream oss;
                oss << "ERROR: Failed to parse fragment range from bam record. Frag begin,end: " << fragBegin << " " << fragEnd << " bamRecord: " << bamRead << "\n";
                BOOST_THROW_EXCEPTION(LogicException(oss.str()));
            }

            const pos_t fragOverlap(std::min((centerPos-fragBegin), (fragEnd-centerPos)));

            if (fragOverlap < minFragSupport) continue;

            SVFragmentEvidence& fragment(evidence.getSample(isTumor)[bamRead.qname()]);
            SVFragmentEvidenceAllele& ref(fragment.ref);

            setReadEvidence(minMapQ, bamRead, fragment.getRead(bamRead.is_first()));

            if (isBp1)
            {
                ref.bp1.isFragmentSupport = true;
            }
            else
            {
                ref.bp2.isFragmentSupport = true;
            }


            if (isDoubleCountSkip) continue;

            /// old tracker:
            if (isBp1)
            {
                sample.ref.bp1SpanReadCount++;
            }
            else
            {
                sample.ref.bp2SpanReadCount++;
            }
        }
    }
}



// make final interpretation of reference support as the minimum breakend support:
static
void
finishSamplePairSupport(
    SVSampleInfo& sample)
{
    sample.ref.spanPairCount = std::min(sample.ref.bp1SpanReadCount, sample.ref.bp2SpanReadCount);
}



void
SVScorer::
getSVRefPairSupport(
    const SVCandidate& sv,
    SVScoreInfo& baseInfo,
    SVEvidence& evidence)
{
    getSVRefPairSupport(sv.bp1, baseInfo, evidence, true);
    getSVRefPairSupport(sv.bp2, baseInfo, evidence, false);

    finishSamplePairSupport(baseInfo.tumor);
    finishSamplePairSupport(baseInfo.normal);
}



void
SVScorer::
getSVPairSupport(
    const SVCandidateSetData& svData,
    const SVCandidate& sv,
    SVScoreInfo& baseInfo,
    SVEvidence& evidence)
{
    const unsigned minMapQ(_readScanner.getMinMapQ());

    // count the read pairs supporting the alternate allele in each sample, using data we already produced during candidate generation:
    //
    const unsigned bamCount(_bamStreams.size());
    for (unsigned bamIndex(0); bamIndex < bamCount; ++bamIndex)
    {
        const bool isTumor(_isAlignmentTumor[bamIndex]);
        SVSampleInfo& sample(isTumor ? baseInfo.tumor : baseInfo.normal);

        const SVCandidateSetReadPairSampleGroup& svDataGroup(svData.getDataGroup(bamIndex));
        BOOST_FOREACH(const SVCandidateSetReadPair& pair, svDataGroup)
        {
            // is this read pair associated with this candidateIndex? (each read pair can be associated with multiple candidates)
            if (0 == std::count(pair.svIndex.begin(),pair.svIndex.end(), sv.candidateIndex)) continue;

            if (! (pair.read1.isSet() || pair.read2.isSet())) continue;

            std::string qname;

            if (pair.read1.isSet())
            {
                qname = pair.read1.bamrec.qname();
            }
            else
            {
                qname = pair.read2.bamrec.qname();
            }

            SVFragmentEvidence& fragment(evidence.getSample(isTumor)[qname]);
            SVFragmentEvidenceAllele& alt(fragment.alt);

            // for all large spanning events -- we don't test for pair support of the two breakends separately -- this could be
            // beneficial if there was an unusually large insertion associated with the event. For now we approximate that
            // these events will mostly not have very large insertions.
            //
            alt.bp1.isFragmentSupport = true;
            alt.bp2.isFragmentSupport = true;

            if (pair.read1.isSet())
            {
                sample.alt.bp1SpanReadCount += 1;
                setReadEvidence(minMapQ, pair.read1.bamrec, fragment.read1);
            }

            if (pair.read2.isSet())
            {
                sample.alt.bp2SpanReadCount += 1;
                setReadEvidence(minMapQ, pair.read2.bamrec, fragment.read2);
            }

            if (pair.read1.isSet() && pair.read2.isSet())
            {
                sample.alt.spanPairCount += 1;
            }

        }
    }

    // count the read pairs supporting the reference allele on each breakend in each sample:
    //
    getSVRefPairSupport(sv, baseInfo, evidence);
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
    //
    // consider 2-locus events first
    // TODO: to add local assembly later
    //
    const bool isSkipSRSearch(
        (! assemblyData.isSpanning) ||
        (sv.isImprecise()) ||
        (isSkipSRSearchDepth));

    if (isSkipSRSearch) return;

    // Get Data on standard read pairs crossing the two breakends,

    // extract SV alignment info for split read evidence
    const SVAlignmentInfo SVAlignInfo(sv, assemblyData);

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
    log_os << "tumor contig SP count: " << baseInfo.tumor.contigSRCount << "\n";
    log_os << "tumor contig SP evidence: " << baseInfo.tumor.contigSREvidence << "\n";
    log_os << "tumor contig SP_mapQ: " << baseInfo.tumor.contigSRMapQ << "\n";
    log_os << "normal contig SP count: " << baseInfo.normal.contigSRCount << "\n";
    log_os << "normal contig SP evidence: " << baseInfo.normal.contigSREvidence << "\n";
    log_os << "normal contig SP_mapQ: " << baseInfo.normal.contigSRMapQ << "\n";

    log_os << "tumor ref SP count: " << baseInfo.tumor.refSRCount << "\n";
    log_os << "tumor ref SP evidence: " << baseInfo.tumor.refSREvidence << "\n";
    log_os << "tumor ref SP_mapQ: " << baseInfo.tumor.refSRMapQ << "\n";
    log_os << "normal ref SP count: " << baseInfo.normal.refSRCount << "\n";
    log_os << "normal ref SP evidence: " << baseInfo.normal.refSREvidence << "\n";
    log_os << "normal ref SP_mapQ: " << baseInfo.normal.refSRMapQ << "\n";
#endif
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

#ifdef DEBUG_SVS
    log_os << SVAlignInfo << "\n";
#endif

    // count the paired read fragments supporting the ref and alt alleles in each sample:
    //
    getSVPairSupport(svData, sv, baseInfo, evidence);

    // count the split reads supporting the ref and alt alleles in each sample
    //
    getSVSplitReadSupport(assemblyData, sv, baseInfo, evidence);
}



/// score diploid germline specific components:
static
void
scoreDiploidSV(
    const CallOptionsDiploid& diploidOpt,
    const SVCandidate& sv,
    const ChromDepthFilterUtil& dFilter,
    SVScoreInfo& baseInfo,
    SVScoreInfoDiploid& diploidInfo)
{
    //
    // compute qualities
    //
    {
        diploidInfo.altScore=60;
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
                baseInfo.filters.insert(diploidOpt.maxDepthFilterLabel);
            }
            else if (baseInfo.bp2MaxDepth > dFilter.maxDepth(sv.bp2.interval.tid))
            {
                baseInfo.filters.insert(diploidOpt.maxDepthFilterLabel);
            }
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
        if (baseInfo.normal.alt.spanPairCount > 1) isNonzeroSomaticQuality=false;

        if (isNonzeroSomaticQuality)
        {
            const bool lowPairSupport(baseInfo.tumor.alt.spanPairCount < 6);
            const bool lowSingleSupport((baseInfo.tumor.alt.bp1SpanReadCount < 14) || (baseInfo.tumor.alt.bp2SpanReadCount < 14));
            const bool highSingleContam((baseInfo.normal.alt.bp1SpanReadCount > 1) || (baseInfo.normal.alt.bp2SpanReadCount > 1));

            /// allow single pair support to rescue an SV only if the evidence looks REALLY good:
            if (lowPairSupport && (lowSingleSupport || highSingleContam))
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
                baseInfo.filters.insert(somaticOpt.maxDepthFilterLabel);
            }
            else if (baseInfo.bp2MaxDepth > dFilter.maxDepth(sv.bp2.interval.tid))
            {
                baseInfo.filters.insert(somaticOpt.maxDepthFilterLabel);
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
    scoreDiploidSV(_diploidOpt, sv, _dFilterDiploid, modelScoreInfo.base, modelScoreInfo.diploid);

    // score components specific to somatic model:
    if (isSomatic)
    {
        scoreSomaticSV(_somaticOpt, sv, _dFilterSomatic, modelScoreInfo.base, modelScoreInfo.somatic);
    }
}




