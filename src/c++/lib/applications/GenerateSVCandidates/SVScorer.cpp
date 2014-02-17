// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Manta
// Copyright (c) 2013-2014 Illumina, Inc.
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
#include "blt_util/math_util.hh"
#include "blt_util/prob_util.hh"
#include "blt_util/qscore.hh"
#include "common/Exceptions.hh"
#include "manta/ReadGroupStatsSet.hh"
#include "manta/SVCandidateUtil.hh"

#include "boost/array.hpp"
#include "boost/foreach.hpp"

#include <algorithm>
#include <iostream>
#include <string>


//#define DEBUG_SCORE
//#define DEBUG_SOMATIC_SCORE

#if defined(DEBUG_SCORE) || defined(DEBUG_SOMATIC_SCORE)
#include "blt_util/log.hh"
#endif




SVScorer::
SVScorer(
    const GSCOptions& opt,
    const bam_header_info& header) :
    _isAlignmentTumor(opt.alignFileOpt.isAlignmentTumor),
    _callOpt(opt.callOpt),
    _callDopt(_callOpt),
    _diploidOpt(opt.diploidOpt),
    _diploidDopt(_diploidOpt),
    _somaticOpt(opt.somaticOpt),
    _somaticDopt(_somaticOpt),
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

        if (is_segment_align_match(ps.type))
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



void
SVScorer::
getBreakendMaxMappedDepthAndMQ0(
    const bool isMaxDepth,
    const double cutoffDepth,
    const SVBreakend& bp,
    unsigned& maxDepth,
    float& MQ0Frac)
{
    /// define a new interval -/+ 50 bases around the center pos
    /// of the breakpoint
    static const pos_t regionSize(50);

    maxDepth=0;
    MQ0Frac=0;

    unsigned totalReads(0);
    unsigned totalMQ0Reads(0);

    const pos_t centerPos(bp.interval.range.center_pos());
    const known_pos_range2 searchRange(std::max((centerPos-regionSize),0), (centerPos+regionSize));

    if (searchRange.size() == 0) return;

    std::vector<unsigned> depth(searchRange.size(),0);

    bool isCutoff(false);
    bool isNormalFound(false);

    const unsigned bamCount(_bamStreams.size());
    for (unsigned bamIndex(0); bamIndex < bamCount; ++bamIndex)
    {
        if (_isAlignmentTumor[bamIndex]) continue;
        isNormalFound=true;

        bam_streamer& bamStream(*_bamStreams[bamIndex]);

        // set bam stream to new search interval:
        bamStream.set_new_region(bp.interval.tid, searchRange.begin_pos(), searchRange.end_pos());

        while (bamStream.next())
        {
            const bam_record& bamRead(*(bamStream.get_record_ptr()));

            // turn filtration down to mapped only to match depth estimate method:
            if (bamRead.is_unmapped()) continue;

            const pos_t refPos(bamRead.pos()-1);
            if (refPos >= searchRange.end_pos()) break;

            addReadToDepthEst(bamRead,searchRange.begin_pos(),depth);

            totalReads++;
            if (0 == bamRead.map_qual()) totalMQ0Reads++;

            if (isMaxDepth)
            {
                const pos_t depthOffset(refPos-searchRange.begin_pos());
                if (depthOffset>=0)
                {
                    if (depth[depthOffset] > cutoffDepth)
                    {
                        isCutoff=true;
                        break;
                    }
                }
            }
        }

        if (isCutoff) break;
    }

    assert(isNormalFound);

    maxDepth = *(std::max_element(depth.begin(),depth.end()));
    if (totalReads>=10)
    {
        MQ0Frac = static_cast<float>(totalMQ0Reads)/static_cast<float>(totalReads);
    }
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
addConservativeSplitReadSupport(
    const SVFragmentEvidence& fragev,
    const bool isRead1,
    SVSampleInfo& sampleBaseInfo)
{
    static const float splitSupportProb(0.999);

    // only consider reads where at least one allele and one breakend is confident
    //
    // ...note this is done in the absence of having a noise state in the model
    //
    if (! fragev.isAnySplitReadSupport(isRead1)) return;

    float altLnLhood =
        std::max(fragev.alt.bp1.getRead(isRead1).splitLnLhood,
                 fragev.alt.bp2.getRead(isRead1).splitLnLhood);

    float refLnLhood =
        std::max(fragev.ref.bp1.getRead(isRead1).splitLnLhood,
                 fragev.ref.bp2.getRead(isRead1).splitLnLhood);

    // convert to normalized prob:
    if (altLnLhood > refLnLhood)
    {
        lnToProb(refLnLhood, altLnLhood);
        if (altLnLhood > splitSupportProb) sampleBaseInfo.alt.confidentSplitReadCount++;
    }
    else
    {
        lnToProb(altLnLhood, refLnLhood);
        if (refLnLhood > splitSupportProb) sampleBaseInfo.ref.confidentSplitReadCount++;
    }
}



static
float
getSpanningPairAlleleLhood(
    const SVFragmentEvidenceAllele& allele)
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

    return fragProb;
}



static
void
addConservativeSpanningPairSupport(
    const SVFragmentEvidence& fragev,
    SVSampleInfo& sampleBaseInfo)
{
    static const float pairSupportProb(0.9);

    if (! fragev.isAnySpanningPairSupport()) return;

    float altLhood(getSpanningPairAlleleLhood(fragev.alt));
    float refLhood(getSpanningPairAlleleLhood(fragev.ref));

    assert(altLhood >= 0);
    assert(refLhood >= 0);
    if ((altLhood <= 0) && (refLhood <= 0))
    {
        using namespace illumina::common;

        std::ostringstream oss;
        oss << "ERROR: Spanning likelihood is zero for all alleles. Fragment: " << fragev << "\n";
        BOOST_THROW_EXCEPTION(LogicException(oss.str()));
    }

    const bool isFullyMapped(fragev.read1.isObservedAnchor() && fragev.read2.isObservedAnchor());

    // convert to normalized prob:
    const float sum(altLhood+refLhood);
    if (altLhood > refLhood)
    {
        if ((altLhood/sum) > pairSupportProb)
        {
            sampleBaseInfo.alt.confidentSemiMappedSpanningPairCount++;
            if (isFullyMapped) sampleBaseInfo.alt.confidentSpanningPairCount++;
        }
    }
    else
    {
        if ((refLhood/sum) > pairSupportProb)
        {
            sampleBaseInfo.ref.confidentSemiMappedSpanningPairCount++;
            if (isFullyMapped) sampleBaseInfo.ref.confidentSpanningPairCount++;
        }
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
        addConservativeSplitReadSupport(fragev,true,sampleBaseInfo);
        addConservativeSplitReadSupport(fragev,false,sampleBaseInfo);

        addConservativeSpanningPairSupport(fragev, sampleBaseInfo);
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
    // at what factor above the maxDepth FILTER criteria do we stop enumerating scoring components?
    static const unsigned cutoffDepthFactor(2);

    const bool isMaxDepth(_dFilterDiploid.isMaxDepthFilter() && _dFilterSomatic.isMaxDepthFilter());
    double bp1CutoffDepth(0);
    double bp2CutoffDepth(0);
    if (isMaxDepth)
    {
        const double bp1MaxMaxDepth(std::max(_dFilterDiploid.maxDepth(sv.bp1.interval.tid), _dFilterSomatic.maxDepth(sv.bp1.interval.tid)));
        const double bp2MaxMaxDepth(std::max(_dFilterDiploid.maxDepth(sv.bp2.interval.tid), _dFilterSomatic.maxDepth(sv.bp2.interval.tid)));

        bp1CutoffDepth = cutoffDepthFactor*bp1MaxMaxDepth;
        bp2CutoffDepth = cutoffDepthFactor*bp2MaxMaxDepth;
    }

    // get breakend center_pos depth estimate:
    getBreakendMaxMappedDepthAndMQ0(isMaxDepth, bp1CutoffDepth, sv.bp1, baseInfo.bp1MaxDepth, baseInfo.bp1MQ0Frac);
    const bool isBp1OverDepth(baseInfo.bp1MaxDepth > bp1CutoffDepth);
    if (! (isMaxDepth && isBp1OverDepth))
    {
        getBreakendMaxMappedDepthAndMQ0(isMaxDepth, bp2CutoffDepth, sv.bp2, baseInfo.bp2MaxDepth, baseInfo.bp2MQ0Frac);
    }
    const bool isBp2OverDepth(baseInfo.bp2MaxDepth > bp2CutoffDepth);

    const bool isOverDepth(isBp1OverDepth || isBp2OverDepth);
    const bool isSkipEvidenceSearch(isMaxDepth && isOverDepth);

    if (! isSkipEvidenceSearch)
    {
        // count the paired-read fragments supporting the ref and alt alleles in each sample:
        //
        getSVPairSupport(svData, assemblyData, sv, evidence);

        // count the split reads supporting the ref and alt alleles in each sample
        //
        getSVSplitReadSupport(assemblyData, sv, baseInfo, evidence);
    }

    // compute allele likelihoods, and any other summary metric shared between all models:
    //
    getSVSupportSummary(evidence, baseInfo);
}



/// record a set of convenient companion values for any probability
///
struct ProbSet
{
    ProbSet(const double initProb) :
        prob(initProb),
        comp(1-prob),
        lnProb(std::log(prob)),
        lnComp(std::log(comp))
    {}

    double prob;
    double comp;
    double lnProb;
    double lnComp;
};



static
void
incrementSpanningPairAlleleLnLhood(
    const ProbSet& chimeraProb,
    const SVFragmentEvidenceAllele& allele,
    const double power,
    double& bpLnLhood)
{
    const float fragProb(getSpanningPairAlleleLhood(allele));
    bpLnLhood += std::log(chimeraProb.comp*fragProb + chimeraProb.prob)*power;
}



static
void
incrementAlleleSplitReadLhood(
    const ProbSet& selfMapProb,
    const ProbSet& otherMapProb,
    const SVFragmentEvidenceAllele& allele,
    const double /*readLnPrior*/,
    const bool isRead1,
    double& refSplitLnLhood,
    bool& isReadEvaluated)
{
    if (! (allele.bp1.getRead(isRead1).isSplitEvaluated &&
           allele.bp2.getRead(isRead1).isSplitEvaluated))
    {
        isReadEvaluated = false;
    }

    const double alignBp1LnLhood(allele.bp1.getRead(isRead1).splitLnLhood);
    const double alignBp2LnLhood(allele.bp2.getRead(isRead1).splitLnLhood);

    const double alignLnLhood(std::max(alignBp1LnLhood,alignBp2LnLhood));

    const double fragLnLhood = log_sum((selfMapProb.lnComp+alignLnLhood), (otherMapProb.lnProb)); //+readLnPrior));
    refSplitLnLhood += fragLnLhood;

#ifdef DEBUG_SCORE
    static const std::string logtag("incrementAlleleSplitReadLhood: ");
    log_os << logtag //<< "readPrior: " << readLnPrior
           << " isRead1?: " << isRead1 << "\n";
    log_os << logtag << "isEval " << isReadEvaluated << "\n";
    log_os << logtag << "alignBp1LnLhood " << alignBp1LnLhood << "\n";
    log_os << logtag << "alignBp2LnLhood " << alignBp2LnLhood << "\n";
    log_os << logtag << "selfMap " << selfMapProb.lnProb << "\n";
    log_os << logtag << "otherMap " << otherMapProb.lnProb << "\n";
    log_os << logtag << "increment " << fragLnLhood << "\n";
    log_os << logtag << "refSplitLnLhood " << refSplitLnLhood << "\n";
#endif

}



static
void
incrementSplitReadLhood(
    const std::string& /*fragLabel*/,
    const SVFragmentEvidence& fragev,
    const ProbSet& refMapProb,
    const ProbSet& altMapProb,
    const bool isPermissive,
    const bool isRead1,
    double& refSplitLnLhood,
    double& altSplitLnLhood,
    bool& isReadEvaluated)
{
    static const double baseLnPrior(std::log(0.25));

#ifdef DEBUG_SCORE
    log_os << __FUNCTION__ << ": pre-support\n";
#endif

    bool isSupported(false);
    if (isPermissive)
    {
        isSupported=fragev.isAnyTier2SplitReadSupport(isRead1);
    }
    else
    {
        isSupported=fragev.isAnySplitReadSupport(isRead1);
    }

    if (! isSupported)
    {
        isReadEvaluated = false;
        return;
    }

#ifdef DEBUG_SCORE
    log_os << __FUNCTION__ << ": post-support\n";
#endif

    const unsigned readSize(fragev.getRead(isRead1).size);
    const double readLnPrior(baseLnPrior*readSize);

#ifdef DEBUG_SCORE
    log_os << __FUNCTION__ << ": starting ref\n";
#endif
    incrementAlleleSplitReadLhood(refMapProb, altMapProb, fragev.ref, readLnPrior, isRead1, refSplitLnLhood, isReadEvaluated);
#ifdef DEBUG_SCORE
    log_os << __FUNCTION__ << ": starting alt\n";
#endif
    incrementAlleleSplitReadLhood(altMapProb, refMapProb, fragev.alt, readLnPrior, isRead1, altSplitLnLhood, isReadEvaluated);
}



struct AlleleLnLhood
{
    AlleleLnLhood() :
        fragPair(0),
        read1Split(0),
        read2Split(0)
    {}

    double fragPair;
    double read1Split;
    double read2Split;
};



static
double
getFragLnLhood(
    const AlleleLnLhood& al,
    const bool isRead1Evaluated,
    const bool isRead2Evaluated)
{
#ifdef DEBUG_SCORE
    log_os << "getFragLnLhood: frag/read1/read2 " << al.fragPair << " " << al.read1Split << " " << al.read2Split << "\n";
    log_os << "getFragLnLhood: isread1/isread2 " << isRead1Evaluated << " " << isRead2Evaluated << "\n";
#endif

    double ret(al.fragPair);

    // limit split read evidence to only one read, b/c it's only possible for one section
    // of the molecule to independently cross the breakend:
    if (isRead1Evaluated)
    {
        if (isRead2Evaluated)
        {
            ret += std::max(al.read1Split, al.read2Split);
        }
        else
        {
            ret += al.read1Split;
        }
    }
    else if (isRead2Evaluated)
    {
        ret += al.read2Split;
    }

    return ret;
}



/// when an sv is treated as 'small', we skip all paired-read evidence and rely on split reads only:
///
/// with further model improvements we can add pairs back into the small variant calls:
///
/// this function returns 1 for a variant which is "fully small" and 0 for a variant which is "fully large",
/// the treatment of paired-read evidence at intermediate values is defined in client scoring code
///
static
float
getSmallSVWeight(
    const SVCandidate& sv,
    const bool isSmallAssembler)
{
    static const int minSmallSize(300);
    static const int maxSmallSize(500);

    assert(maxSmallSize >= minSmallSize);
    static const int smallRange(maxSmallSize-minSmallSize);

    if (! isSmallAssembler) return 0.;
    const int svSize(std::abs(sv.bp2.interval.range.center_pos() - sv.bp1.interval.range.center_pos()));
    return std::min(1.f,std::max(0.f,static_cast<float>(maxSmallSize-svSize)/static_cast<float>(smallRange)));
}



static
float
largeNoiseSVPriorWeight(
    const SVCandidate& sv)
{
    static const int smallSize(5000);
    static const int largeSize(10000);

    const int svSize(std::abs(sv.bp2.interval.range.center_pos() - sv.bp1.interval.range.center_pos()));
    return std::min(static_cast<float>(std::max(svSize-smallSize,0))/static_cast<float>(largeSize),1.f);
}



/// return true if any evidence exists for fragment:
///
/// \param semiMappedPower multiply out semi-mapped reads (in log space) by this value
///
static
bool
getRefAltFromFrag(
    const float smallSVWeight,
    const double semiMappedPower,
    const ProbSet& chimeraProb,
    const ProbSet& refSplitMapProb,
    const ProbSet& altSplitMapProb,
    const bool isPermissive,
    const std::string& fragLabel,
    const SVFragmentEvidence& fragev,
    AlleleLnLhood& refLnLhoodSet,
    AlleleLnLhood& altLnLhoodSet,
    bool& isRead1Evaluated,
    bool& isRead2Evaluated)
{
#ifdef DEBUG_SCORE
    log_os << __FUNCTION__ << ": qname: " << fragLabel << " fragev: " << fragev << "\n";
#endif

    /// TODO: add read pairs with one shadow read to the alt read pool

    /// high-quality spanning support relies on read1 and read2 mapping well:
    bool isFragEvaluated(false);

    const bool isPairUsable = ((fragev.read1.isScanned && fragev.read2.isScanned) &&
                               (fragev.read1.isAnchored || fragev.read2.isAnchored));

    if (isPairUsable)
    {
        /// only add to the likelihood if the fragment "supports" at least one allele:
        if ( fragev.isAnySpanningPairSupport() )
        {
            {
                const bool isSemiMapped(! (fragev.read1.isAnchored && fragev.read2.isAnchored));
                double spanPower( isSemiMapped ?  semiMappedPower : 1.);
                spanPower *= (1.-smallSVWeight);

                incrementSpanningPairAlleleLnLhood(chimeraProb, fragev.ref, spanPower, refLnLhoodSet.fragPair);
                incrementSpanningPairAlleleLnLhood(chimeraProb, fragev.alt, spanPower, altLnLhoodSet.fragPair);
                isFragEvaluated=true;
            }
        }
    }

    /// split support is less dependent on mapping quality of the individual read, because
    /// we're potentially relying on shadow reads recovered from the unmapped state
    isRead1Evaluated = true;
    isRead2Evaluated = true;
#ifdef DEBUG_SCORE
    log_os << __FUNCTION__ << ": starting read1 split\n";
#endif
    incrementSplitReadLhood(fragLabel, fragev, refSplitMapProb, altSplitMapProb, isPermissive, true,  refLnLhoodSet.read1Split, altLnLhoodSet.read1Split, isRead1Evaluated);
#ifdef DEBUG_SCORE
    log_os << __FUNCTION__ << ": starting read2 split\n";
#endif
    incrementSplitReadLhood(fragLabel, fragev, refSplitMapProb, altSplitMapProb, isPermissive, false, refLnLhoodSet.read2Split, altLnLhoodSet.read2Split, isRead2Evaluated);

#ifdef DEBUG_SCORE
    log_os << __FUNCTION__ << ": iseval frag/read1/read2: " << isFragEvaluated << " " << isRead1Evaluated << " " << isRead1Evaluated << "\n";
#endif
    return (isFragEvaluated || isRead1Evaluated || isRead2Evaluated);
}



/// score diploid germline specific components:
static
void
addDiploidLoglhood(
    const float smallSVWeight,
    const SVEvidence::evidenceTrack_t& sampleEvidence,
    boost::array<double,DIPLOID_GT::SIZE>& loglhood)
{
    BOOST_FOREACH(const SVEvidence::evidenceTrack_t::value_type& val, sampleEvidence)
    {
        const std::string& fragLabel(val.first);
        const SVFragmentEvidence& fragev(val.second);

        AlleleLnLhood refLnLhoodSet, altLnLhoodSet;
        bool isRead1Evaluated(true);
        bool isRead2Evaluated(true);

        /// TODO: set this from graph data:
        ///
        /// put some more thought into this -- is this P (spurious | any old read) or P( spurious | chimera ) ??
        /// it seems like it should be the latter in the usages that really matter.
        ///
        static const ProbSet chimeraProb(1e-3);

        /// use a constant mapping prob for now just to get the zero-th order concept into the model
        /// that "reads are mismapped at a non-trivial rate"
        /// TODO: experiment with per-read mapq values
        ///
        static const ProbSet refSplitMapProb(1e-6);
        static const ProbSet altSplitMapProb(1e-4);

        /// don't use semi-mapped reads for germline calling:
        static const double semiMappedPower(0.);

        static const bool isPermissive(false);

        if (! getRefAltFromFrag(smallSVWeight, semiMappedPower, chimeraProb, refSplitMapProb, altSplitMapProb, isPermissive, fragLabel, fragev,
                                refLnLhoodSet, altLnLhoodSet, isRead1Evaluated, isRead2Evaluated))
        {
            // continue if this fragment was not evaluated for pair or split support for either allele:
            continue;
        }

        for (unsigned gt(0); gt<DIPLOID_GT::SIZE; ++gt)
        {
            using namespace DIPLOID_GT;

#ifdef DEBUG_SCORE
            log_os << __FUNCTION__ << ": starting gt: " << gt << " " << label(gt) << "\n";
#endif

            const index_t gtid(static_cast<const index_t>(gt));
            const double refLnFragLhood(getFragLnLhood(refLnLhoodSet, isRead1Evaluated, isRead2Evaluated));
#ifdef DEBUG_SCORE
            log_os << __FUNCTION__ << ": refLnFragLhood: " << refLnFragLhood << "\n";
#endif
            const double altLnFragLhood(getFragLnLhood(altLnLhoodSet, isRead1Evaluated, isRead2Evaluated));
#ifdef DEBUG_SCORE
            log_os << __FUNCTION__ << ": altLnFragLhood: " << altLnFragLhood << "\n";
#endif
            const double refLnLhood(refLnFragLhood + altLnCompFraction(gtid));
            const double altLnLhood(altLnFragLhood + altLnFraction(gtid));
            loglhood[gt] += log_sum(refLnLhood, altLnLhood);

#ifdef DEBUG_SCORE
            log_os << __FUNCTION__ << ": gt/fragref/ref/fragalt/alt: "
                   << label(gt)
                   << " " << refLnFragLhood
                   << " " << refLnLhood
                   << " " << altLnFragLhood
                   << " " << altLnLhood
                   << "\n";
#endif
        }
    }
}




/// score diploid germline specific components:
static
void
scoreDiploidSV(
    const CallOptionsDiploid& diploidOpt,
    const CallOptionsDiploidDeriv& diploidDopt,
    const SVCandidate& sv,
    const float smallSVWeight,
    const ChromDepthFilterUtil& dFilter,
    const SVEvidence& evidence,
    SVScoreInfo& baseInfo,
    SVScoreInfoDiploid& diploidInfo)
{
    //
    // compute qualities
    //
    {
        boost::array<double,DIPLOID_GT::SIZE> loglhood;
        std::fill(loglhood.begin(),loglhood.end(),0);
        addDiploidLoglhood(smallSVWeight, evidence.normal, loglhood);

        boost::array<double,DIPLOID_GT::SIZE> pprob;
        for (unsigned gt(0); gt<DIPLOID_GT::SIZE; ++gt)
        {
            // TODO: fix from xchen to use log prior instead of prior
            // swap this in once there's time to make any compensatory parameter adjustments
            //pprob[gt] = loglhood[gt] + diploidDopt.logPrior[gt];

            // this version is wrong, it uses the prior without the log, see above:
            pprob[gt] = loglhood[gt] + diploidDopt.prior[gt];
        }

        unsigned maxGt(0);
        normalize_ln_distro(pprob.begin(), pprob.end(), maxGt);

#ifdef DEBUG_SCORE
        for (unsigned gt(0); gt<DIPLOID_GT::SIZE; ++gt)
        {
            log_os << __FUNCTION__ << ": gt/lhood/prior/pprob: "
                   << DIPLOID_GT::label(gt)
                   << " " << loglhood[gt]
                   << " " << diploidDopt.prior[gt]
                   << " " << pprob[gt]
                   << "\n";
        }
#endif

        diploidInfo.gt=static_cast<DIPLOID_GT::index_t>(maxGt);
        diploidInfo.altScore=error_prob_to_qphred(pprob[DIPLOID_GT::REF]);
        diploidInfo.gtScore=error_prob_to_qphred(prob_comp(pprob.begin(),pprob.end(), diploidInfo.gt));
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

        if (diploidInfo.gtScore < diploidOpt.minPassGTScore)
        {
            diploidInfo.filters.insert(diploidOpt.minGTFilterLabel);
        }

        const bool isMQ0FilterSize(isSVBelowMinSize(sv,1000));
        if (isMQ0FilterSize)
        {
            if ((baseInfo.bp1MQ0Frac > diploidOpt.maxMQ0Frac) ||
                (baseInfo.bp2MQ0Frac > diploidOpt.maxMQ0Frac))
            {
                diploidInfo.filters.insert(diploidOpt.maxMQ0FracLabel);
            }
        }
    }
}



static
unsigned
getSpanningPairCount(
    const SVSampleAlleleInfo& allele,
    const float smallSVWeight,
    const bool isPermissive)
{
    if (isPermissive) return (1.-smallSVWeight)*allele.confidentSemiMappedSpanningPairCount;
    else              return (1.-smallSVWeight)*allele.confidentSpanningPairCount;
}



static
unsigned
getSupportCount(
    const SVSampleAlleleInfo& allele,
    const float smallSVWeight,
    const bool isPermissive)
{
    return allele.confidentSplitReadCount + getSpanningPairCount(allele, smallSVWeight, isPermissive);
}



static
double
estimateSomaticMutationFreq(
    const SVScoreInfo& baseInfo,
    const float smallSVWeight,
    const bool isPermissive)
{
    const unsigned altCounts = getSupportCount(baseInfo.tumor.alt, smallSVWeight, isPermissive);
    const unsigned refCounts = getSupportCount(baseInfo.tumor.ref, smallSVWeight, isPermissive);
    if ((altCounts + refCounts) == 0) return 0;
    return static_cast<double>(altCounts) / static_cast<double>(altCounts + refCounts);
}



static
double
estimateNoiseMutationFreq(
    const SVScoreInfo& baseInfo,
    const float smallSVWeight,
    const bool isPermissive)
{
    const unsigned normalAltCounts = getSupportCount(baseInfo.normal.alt, smallSVWeight, isPermissive);
    const unsigned normalRefCounts = getSupportCount(baseInfo.normal.ref, smallSVWeight, isPermissive);
    const unsigned tumorAltCounts = getSupportCount(baseInfo.tumor.alt, smallSVWeight, isPermissive);
    const unsigned tumorRefCounts = getSupportCount(baseInfo.tumor.ref, smallSVWeight, isPermissive);

    const unsigned altCounts(normalAltCounts + tumorAltCounts);
    const unsigned refCounts(normalRefCounts + tumorRefCounts);

    if ((altCounts + refCounts) == 0) return 0;
    return static_cast<double>(altCounts) / static_cast<double>(altCounts + refCounts);
}



static
void
computeSomaticSampleLoghood(
    const float smallSVWeight,
    const SVEvidence::evidenceTrack_t& evidenceTrack,
    const double somaticMutationFreq,
    const double noiseMutationFreq,
    const double isPermissive,
    boost::array<double,SOMATIC_GT::SIZE>& loglhood)
{
    /// TODO: find a better way to set this number from training data:
    static const ProbSet chimeraProbDefault(1e-4);
    static const ProbSet chimeraProbPermissive(1e-5);
    const ProbSet& chimeraProb( isPermissive ? chimeraProbPermissive : chimeraProbDefault );

    /// use a constant mapping prob for now just to get the zero-th order concept into the model
    /// that "reads are mismapped at a non-trivial rate"
    /// TODO: experiment with per-read mapq values
    ///
    static const ProbSet refSplitMapProb(1e-6);

    static const ProbSet altSplitMapProbDefault(1e-4);
    static const ProbSet altSplitMapProbPermissive(1e-6);
    const ProbSet& altSplitMapProb( isPermissive ? altSplitMapProbPermissive : altSplitMapProbDefault );

    // semi-mapped reads make a partial contribution in tier1, and a full contribution in tier2:
    const double semiMappedPower( isPermissive ? 1. : 0. );

    BOOST_FOREACH(const SVEvidence::evidenceTrack_t::value_type& val, evidenceTrack)
    {
        const std::string& fragLabel(val.first);
        const SVFragmentEvidence& fragev(val.second);

        AlleleLnLhood refLnLhoodSet, altLnLhoodSet;
        bool isRead1Evaluated(true);
        bool isRead2Evaluated(true);

        if (! getRefAltFromFrag(smallSVWeight, semiMappedPower, chimeraProb, refSplitMapProb, altSplitMapProb, isPermissive, fragLabel, fragev,
                                refLnLhoodSet, altLnLhoodSet, isRead1Evaluated, isRead2Evaluated))
        {
            // continue if this fragment was not evaluated for pair or split support for either allele:
            continue;
        }

        for (unsigned gt(0); gt<SOMATIC_GT::SIZE; ++gt)
        {
            using namespace SOMATIC_GT;

#ifdef DEBUG_SCORE
            log_os << __FUNCTION__ << ": starting gt: " << gt << " " << label(gt) << "\n";
#endif

            const index_t gtid(static_cast<const index_t>(gt));

            const double refLnFragLhood(getFragLnLhood(refLnLhoodSet, isRead1Evaluated, isRead2Evaluated));
#ifdef DEBUG_SCORE
            log_os << __FUNCTION__ << ": refLnFragLhood: " << refLnFragLhood << "\n";
#endif
            const double altLnFragLhood(getFragLnLhood(altLnLhoodSet, isRead1Evaluated, isRead2Evaluated));
#ifdef DEBUG_SCORE
            log_os << __FUNCTION__ << ": altLnFragLhood: " << altLnFragLhood << "\n";
#endif

            // update likelihood with Pr[allele | G]
            const double refLnLhood = refLnFragLhood + altLnCompFraction(gtid, somaticMutationFreq, noiseMutationFreq);
            const double altLnLhood = altLnFragLhood + altLnFraction(gtid, somaticMutationFreq, noiseMutationFreq);

            loglhood[gt] += log_sum(refLnLhood, altLnLhood);
        }
    }
}



/// score somatic specific components:
static
void
scoreSomaticSV(
    const CallOptionsSomatic& somaticOpt,
    const CallOptionsSomaticDeriv& somaticDopt,
    const SVCandidate& sv,
    const float smallSVWeight,
    const ChromDepthFilterUtil& dFilter,
    const SVEvidence& evidence,
    const SVScoreInfo& baseInfo,
    SVScoreInfoSomatic& somaticInfo)
{
    //
    // compute somatic score
    //

    // somatic score is computed at a high stringency date tier (1) and low stringency tier (2), the min value is
    // kept as the final reported quality:
    static const unsigned tierCount(2);
    double tierScore[tierCount] = { 0. , 0. };

    const float largeNoiseWeight(largeNoiseSVPriorWeight(sv));

    for (unsigned tierIndex(0); tierIndex<tierCount; ++tierIndex)
    {
        const bool isPermissive(tierIndex != 0);

        boost::array<double,SOMATIC_GT::SIZE> normalSomaticLhood;
        boost::array<double,SOMATIC_GT::SIZE> tumorSomaticLhood;
        std::fill(normalSomaticLhood.begin(),normalSomaticLhood.end(),0);
        std::fill(tumorSomaticLhood.begin(),tumorSomaticLhood.end(),0);

        // estimate the somatic mutation rate using alternate allele freq from the tumor sample
        const double somaticMutationFreq = estimateSomaticMutationFreq(baseInfo, smallSVWeight, isPermissive);

        // estimate the noise mutation rate using alternate allele freq from the tumor and normal samples
        const double noiseMutationFreq = estimateNoiseMutationFreq(baseInfo, smallSVWeight, isPermissive);

#ifdef DEBUG_SOMATIC_SCORE
        log_os << __FUNCTION__ << ": somaticMutationFrequency: " << somaticMutationFreq << "\n";
        log_os << __FUNCTION__ << ": noiseMutationFrequency: " << noiseMutationFreq << "\n";
        log_os << __FUNCTION__ << ": largeNoiseWeight: " << largeNoiseWeight << "\n";
#endif

        // compute likelihood for the fragments from the tumor sample
        computeSomaticSampleLoghood(smallSVWeight, evidence.tumor, somaticMutationFreq, noiseMutationFreq, isPermissive, tumorSomaticLhood);

        // compute likelihood for the fragments from the normal sample
        computeSomaticSampleLoghood(smallSVWeight, evidence.normal, 0, noiseMutationFreq, isPermissive, normalSomaticLhood);

        boost::array<double,SOMATIC_GT::SIZE> somaticPprob;
        for (unsigned gt(0); gt<SOMATIC_GT::SIZE; ++gt)
        {
            somaticPprob[gt] = tumorSomaticLhood[gt] + normalSomaticLhood[gt] + somaticDopt.logPrior(gt,largeNoiseWeight);
        }

        {
            unsigned maxGt(0);
            normalize_ln_distro(somaticPprob.begin(), somaticPprob.end(), maxGt);
        }

        // independently estimate diploid genotype:
        boost::array<double,DIPLOID_GT::SIZE> normalLhood;
        std::fill(normalLhood.begin(),normalLhood.end(),0);
        addDiploidLoglhood(smallSVWeight, evidence.normal, normalLhood);

        boost::array<double,DIPLOID_GT::SIZE> normalPprob;
        for (unsigned gt(0); gt<DIPLOID_GT::SIZE; ++gt)
        {
            normalPprob[gt] = normalLhood[gt]; // uniform prior for now....
        }

        {
            unsigned maxGt(0);
            normalize_ln_distro(normalPprob.begin(), normalPprob.end(), maxGt);
        }

#ifdef DEBUG_SOMATIC_SCORE
        for (unsigned gt(0); gt<SOMATIC_GT::SIZE; ++gt)
        {
            log_os << __FUNCTION__ << ": somatic gt/tumor_lhood/normal_lhood/prior/pprob: "
                   << SOMATIC_GT::label(gt)
                   << " " << tumorSomaticLhood[gt]
                   << " " << normalSomaticLhood[gt]
                   << " " << somaticDopt.logPrior(gt,largeNoiseWeight)
                   << " " << somaticPprob[gt]
                   << "\n";
        }

        for (unsigned gt(0); gt<DIPLOID_GT::SIZE; ++gt)
        {
            log_os << __FUNCTION__ << ": diploid gt/lhood/pprob: "
                   << DIPLOID_GT::label(gt)
                   << " " << normalLhood[gt]
                   << " " << normalPprob[gt]
                   << "\n";
        }
#endif

        const double nonsomaticProb(prob_comp(somaticPprob.begin(), somaticPprob.end(), SOMATIC_GT::SOM));
        const double nonrefProb(prob_comp(normalPprob.begin(), normalPprob.end(), DIPLOID_GT::REF));

        // not (somatic AND normal ref):
        // (1-(1-a)(1-b)) -> a+b-(ab)
        const double nonsomatic_ref_prob(nonsomaticProb+nonrefProb-(nonsomaticProb*nonrefProb));

        tierScore[tierIndex]=error_prob_to_qphred(nonsomatic_ref_prob);

#ifdef DEBUG_SOMATIC_SCORE
        log_os << __FUNCTION__ << ": tier: " << tierIndex << " somatic score: " << tierScore[tierIndex] << "\n";
#endif

        // don't bother with tier2 if tier1 is too low:
        if (tierScore[tierIndex] <= 0.) break;
    }

    somaticInfo.somaticScore=std::min(tierScore[0],tierScore[1]);

    somaticInfo.somaticScoreTier = 0;
    if (tierScore[1] > tierScore[0])
    {
        somaticInfo.somaticScoreTier = 1;
    }

    // original threshold-based score logic:
#ifdef COUNT_SCORE
    {
        // don't use pair evidence for small variants:
        bool isSmall(false);
        if (isSmallSV)
        {
            const int svSize(std::abs(sv.bp2.interval.range.center_pos() - sv.bp1.interval.range.center_pos()));
            isSmall=(svSize<minPairVariantSize);
        }

#ifdef DEBUG_SCORE
        log_os << __FUNCTION__ << ": isSmall: " << isSmall << '\n'
               << __FUNCTION__ << ": normal: " << baseInfo.normal << '\n'
               << __FUNCTION__ << ": tumor: " << baseInfo.tumor << '\n';
#endif

        bool isNonzeroSomaticQuality(true);

        /// first check for substantial support in the normal:
        if (! isSmall)
        {
            if (baseInfo.normal.alt.confidentSpanningPairCount > 1) isNonzeroSomaticQuality=false;
        }
        if ((baseInfo.normal.alt.confidentSplitReadCount > 0) &&
            ((baseInfo.normal.alt.confidentSplitReadCount > 1) || (baseInfo.tumor.alt.confidentSplitReadCount < 10))) isNonzeroSomaticQuality=false;

        if (isNonzeroSomaticQuality)
        {
            const bool lowPairSupport(baseInfo.tumor.alt.confidentSpanningPairCount < 6);
            const bool lowSplitSupport(baseInfo.tumor.alt.confidentSplitReadCount < 6);
            const bool lowSingleSupport(baseInfo.tumor.alt.confidentSemiMappedSpanningPairCount < 14);
            const bool highSingleContam(baseInfo.normal.alt.confidentSemiMappedSpanningPairCount > 1);

            if (isSmall)
            {
                if (lowSplitSupport) isNonzeroSomaticQuality=false;
            }
            else
            {
                /// allow single pair support to rescue an SV only if the evidence looks REALLY good:
                if ((lowPairSupport && lowSplitSupport) && (lowSingleSupport || highSingleContam))
                    isNonzeroSomaticQuality=false;
            }
        }

        if ((! isSmall) && isNonzeroSomaticQuality)
        {
            if (baseInfo.normal.alt.confidentSpanningPairCount)
            {
                const double ratio(static_cast<double>(baseInfo.tumor.alt.confidentSpanningPairCount)/static_cast<double>(baseInfo.normal.alt.confidentSpanningPairCount));
                if (ratio<9)
                {
                    isNonzeroSomaticQuality=false;
                }
            }
            if (baseInfo.normal.alt.confidentSemiMappedSpanningPairCount)
            {
                const double ratio(static_cast<double>(baseInfo.tumor.alt.confidentSemiMappedSpanningPairCount)/static_cast<double>(baseInfo.normal.alt.confidentSemiMappedSpanningPairCount));
                if (ratio<9)
                {
                    isNonzeroSomaticQuality=false;
                }
            }
        }

        {
            // there needs to be some ref support in the normal as well:
            const bool normRefPairSupport(baseInfo.normal.ref.confidentSpanningPairCount > 6);
            const bool normRefSplitSupport(baseInfo.normal.ref.confidentSplitReadCount > 6);

            if (! (((! isSmall) && normRefPairSupport) || normRefSplitSupport)) isNonzeroSomaticQuality=false;
        }

        if (isNonzeroSomaticQuality) somaticInfo.somaticScore=60;
    }
#endif

    //
    // apply filters
    //
    if (somaticInfo.somaticScore >= somaticOpt.minOutputSomaticScore)
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

        if (somaticInfo.somaticScore < somaticOpt.minPassSomaticScore)
        {
            somaticInfo.filters.insert(somaticOpt.minSomaticScoreLabel);
        }

        const bool isMQ0FilterSize(isSVBelowMinSize(sv,1000));
        if (isMQ0FilterSize)
        {
            if ((baseInfo.bp1MQ0Frac > somaticOpt.maxMQ0Frac) ||
                (baseInfo.bp2MQ0Frac > somaticOpt.maxMQ0Frac))
            {
                somaticInfo.filters.insert(somaticOpt.maxMQ0FracLabel);
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
    const bool isSmallAssembler(! assemblyData.isSpanning);
    const float smallSVWeight(getSmallSVWeight(sv,isSmallAssembler));
    scoreDiploidSV(_diploidOpt, _diploidDopt, sv, smallSVWeight, _dFilterDiploid, evidence, modelScoreInfo.base, modelScoreInfo.diploid);

    // score components specific to somatic model:
    if (isSomatic)
    {
        scoreSomaticSV(_somaticOpt, _somaticDopt, sv, smallSVWeight, _dFilterSomatic, evidence, modelScoreInfo.base, modelScoreInfo.somatic);
    }
}

