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
#include "blt_util/math_util.hh"
#include "blt_util/prob_util.hh"
#include "blt_util/qscore.hh"
#include "common/Exceptions.hh"
#include "manta/ReadGroupStatsSet.hh"

#include "boost/foreach.hpp"

#include <algorithm>
#include <iostream>
#include <string>


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
    _callOpt(opt.callOpt),
    _callDopt(_callOpt),
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
    if (! fragev.isAnySplitSupportForRead(isRead1)) return;

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

    if (! fragev.isAnyPairSupport()) return;

    /// high-quality spanning support relies on read1 and read2 mapping well:
    if (! (fragev.read1.isObservedAnchor() && fragev.read2.isObservedAnchor())) return;

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

    // convert to normalized prob:
    const float sum(altLhood+refLhood);
    if (altLhood > refLhood)
    {
        if ((altLhood/sum) > pairSupportProb) sampleBaseInfo.alt.confidentSpanningPairCount++;
    }
    else
    {
        if ((refLhood/sum) > pairSupportProb) sampleBaseInfo.ref.confidentSpanningPairCount++;
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
    double& bpLnLhood)
{
    const float fragProb(getSpanningPairAlleleLhood(allele));
    bpLnLhood += std::log(chimeraProb.comp*fragProb + chimeraProb.prob);
}



static
void
incrementAlleleSplitReadLhood(
    const SVFragmentEvidenceAllele& allele,
    const double readLnPrior,
    const bool isRead1,
    double& refSplitLnLhood)
{
    /// use a constant mapping prob for now just to get the zero-th order concept into the model
    /// that "reads are mismapped at a non-trivial rate"
    /// TODO: experiment with per-read mapq values here
    ///
    static const double mapProb(1e-6);
    static const double mapComp(1.-mapProb);
    static const double mapLnProb(std::log(mapProb));
    static const double mapLnComp(std::log(mapComp));

    const double alignBp1LnLhood(allele.bp1.getRead(isRead1).splitLnLhood);
    const double alignBp2LnLhood(allele.bp2.getRead(isRead1).splitLnLhood);
    const double alignLnLhood(std::max(alignBp1LnLhood,alignBp2LnLhood));

    refSplitLnLhood += log_sum((mapLnComp+alignLnLhood), (mapLnProb+readLnPrior));

#ifdef DEBUG_SCORE
    static const std::string logtag("incrementAlleleSplitReadLhood: ");
    log_os << logtag << "readPrior: " << readLnPrior << " isRead1?: " << isRead1 << "\n";
    log_os << logtag << "alignBp1LnLhood " << alignBp1LnLhood << "\n";
    log_os << logtag << "alignBp2LnLhood " << alignBp2LnLhood << "\n";
    log_os << logtag << "increment " << log_sum((mapLnComp+alignLnLhood), (mapLnProb+readLnPrior)) << "\n";
    log_os << logtag << "refSplitLnLhood " << refSplitLnLhood << "\n";
#endif

}



static
void
incrementSplitReadLhood(
    const SVFragmentEvidence& fragev,
    const bool isRead1,
    double& refSplitLnLhood,
    double& altSplitLnLhood)
{
    static const double baseLnPrior(std::log(0.25));

#ifdef DEBUG_SCORE
    static const std::string logtag("incrementSplitReadLhood: ");
    log_os << logtag << "pre-support\n";
#endif

    if (! fragev.isAnySplitSupportForRead(isRead1)) return;

#ifdef DEBUG_SCORE
    log_os << logtag << "post-support\n";
#endif

    const unsigned readSize(fragev.getRead(isRead1).size);
    const double readLnPrior(baseLnPrior*readSize);

#ifdef DEBUG_SCORE
    log_os << logtag << "starting ref\n";
#endif
    incrementAlleleSplitReadLhood(fragev.ref, readLnPrior, isRead1, refSplitLnLhood);
#ifdef DEBUG_SCORE
    log_os << logtag << "starting alt\n";
#endif
    incrementAlleleSplitReadLhood(fragev.alt, readLnPrior, isRead1, altSplitLnLhood);
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
    const AlleleLnLhood& al)
{
#ifdef DEBUG_SCORE
    log_os << "getFragLnLhood: frag/read1/read2 " << al.fragPair << " " << al.read1Split << " " << al.read2Split << "\n";
#endif

    return (al.fragPair + al.read1Split + al.read2Split);
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
    static const ProbSet chimeraProb(1e-3);

    //
    // compute qualities
    //
    {
        double loglhood[DIPLOID_GT::SIZE];
        for (unsigned gt(0); gt<DIPLOID_GT::SIZE; ++gt)
        {
            loglhood[gt] = 0.;
        }

        BOOST_FOREACH(const SVEvidence::evidenceTrack_t::value_type& val, evidence.normal)
        {
            const SVFragmentEvidence& fragev(val.second);

            AlleleLnLhood refLnLhoodSet, altLnLhoodSet;

#ifdef DEBUG_SCORE
            log_os << logtag << "qname: " << val.first << " fragev: " << fragev << "\n";
#endif

            /// TODO: add read pairs with one shadow read to the alt read pool

            /// high-quality spanning support relies on read1 and read2 mapping well:
            if ( fragev.read1.isObservedAnchor() && fragev.read2.isObservedAnchor())
            {
                /// only add to the likelihood if the fragment "supports" at least one allele:
                if ( fragev.isAnyPairSupport() )
                {
                    incrementSpanningPairAlleleLnLhood(chimeraProb, fragev.ref, refLnLhoodSet.fragPair);
                    incrementSpanningPairAlleleLnLhood(chimeraProb, fragev.alt, altLnLhoodSet.fragPair);
                }
            }

            /// split support is less dependent on mapping quality of the individual read, because
            /// we're potentially relying on shadow reads recovered from the unmapped state
#ifdef DEBUG_SCORE
            log_os << logtag << "starting read1 split\n";
#endif
            incrementSplitReadLhood(fragev, true,  refLnLhoodSet.read1Split, altLnLhoodSet.read1Split);
#ifdef DEBUG_SCORE
            log_os << logtag << "starting read2 split\n";
#endif
            incrementSplitReadLhood(fragev, false, refLnLhoodSet.read2Split, altLnLhoodSet.read2Split);

            for (unsigned gt(0); gt<DIPLOID_GT::SIZE; ++gt)
            {
                using namespace DIPLOID_GT;

#ifdef DEBUG_SCORE
                log_os << logtag << "starting gt: " << gt << " " << label(gt) << "\n";
#endif

                const index_t gtid(static_cast<const index_t>(gt));
                const double refLnFragLhood(getFragLnLhood(refLnLhoodSet));
#ifdef DEBUG_SCORE
                log_os << logtag << "refLnFragLhood: " << refLnFragLhood << "\n";
#endif
                const double altLnFragLhood(getFragLnLhood(altLnLhoodSet));
#ifdef DEBUG_SCORE
                log_os << logtag << "altLnFragLhood: " << altLnFragLhood << "\n";
#endif
                const double refLnLhood(refLnFragLhood + altLnCompFraction(gtid));
                const double altLnLhood(altLnFragLhood + altLnFraction(gtid));
                loglhood[gt] += log_sum(refLnLhood, altLnLhood);

#ifdef DEBUG_SCORE
                log_os << logtag << "gt/fragref/ref/fragalt/alt: "
                       << label(gt)
                       << " " << refLnFragLhood
                       << " " << refLnLhood
                       << " " << altLnFragLhood
                       << " " << altLnLhood
                       << "\n";
#endif
            }
        }

        double pprob[DIPLOID_GT::SIZE];
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
        if (baseInfo.normal.alt.confidentSpanningPairCount > 1) isNonzeroSomaticQuality=false;
        if (baseInfo.normal.alt.confidentSplitReadCount > 5) isNonzeroSomaticQuality=false;

        if (isNonzeroSomaticQuality)
        {
            const bool lowPairSupport(baseInfo.tumor.alt.confidentSpanningPairCount < 6);
            const bool lowSplitSupport(baseInfo.tumor.alt.confidentSplitReadCount < 6);
            const bool lowSingleSupport((baseInfo.tumor.alt.bp1SpanReadCount < 14) || (baseInfo.tumor.alt.bp2SpanReadCount < 14));
            const bool highSingleContam((baseInfo.normal.alt.bp1SpanReadCount > 1) || (baseInfo.normal.alt.bp2SpanReadCount > 1));

            /// allow single pair support to rescue an SV only if the evidence looks REALLY good:
            if ((lowPairSupport && lowSplitSupport) && (lowSingleSupport || highSingleContam))
                isNonzeroSomaticQuality=false;
        }

        if (isNonzeroSomaticQuality)
        {
            if (baseInfo.normal.alt.confidentSpanningPairCount)
            {
                const double ratio(static_cast<double>(baseInfo.tumor.alt.confidentSpanningPairCount)/static_cast<double>(baseInfo.normal.alt.confidentSpanningPairCount));
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
            const bool normRefPairSupport(baseInfo.normal.ref.confidentSpanningPairCount > 6);
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

