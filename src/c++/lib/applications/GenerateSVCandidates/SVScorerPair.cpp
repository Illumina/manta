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
#include "blt_util/bam_record_util.hh"
#include "common/Exceptions.hh"
#include "manta/SVCandidateUtil.hh"

#include "boost/foreach.hpp"

#include <iostream>
#include <sstream>
#include <string>

/// standard debug output for this file:
//#define DEBUG_PAIR

/// ridiculous debug output for this file:
//#define DEBUG_MEGAPAIR

#ifdef DEBUG_PAIR
#include "blt_util/log.hh"
#endif





static
void
setAlleleFrag(
    const SizeDistribution& fragDistro,
    const int size,
    SVFragmentEvidenceAlleleBreakend& bp)
{
    float fragProb(fragDistro.cdf(size));
    fragProb = std::min(fragProb, (1-fragProb));
#ifdef DEBUG_MEGAPAIR
    log_os << __FUNCTION__ << ": fraglen,prob " << size << " " << fragProb << "\n";
#endif

    bp.isFragmentSupport = true;
    bp.fragLengthProb = fragProb;
}



/// search for all read pairs supporting the alternate allele
///
/// SV types are restricted to be simple insert/delete events and precise (ie. they are all outputs of the small-assembler)
///
/// TODO: adjust this function to understand multiple indels on one haplotype:
///
void
SVScorer::
getSimpleSVAltPairSupport(
    const PairOptions& pairOpt,
    const SVCandidate& sv,
    const bool isBp1,
    SVScoreInfo& baseInfo,
    SVEvidence& evidence)
{
    using namespace illumina::common;

    assert(sv.bp1.interval.tid == sv.bp2.interval.tid);
    assert(getSVType(sv) == SV_TYPE::INDEL);

    /// In case of breakend microhomology approximate the breakend as a point event at the center of the possible range:
    const pos_t centerPos1(sv.bp1.interval.range.center_pos());
    const pos_t centerPos2(sv.bp2.interval.range.center_pos());
    if (centerPos2 <= centerPos1)
    {
        std::ostringstream oss;
        oss << "ERROR: Unexpected breakend orientation in pair support routine for sv: " << sv << "\n";
        BOOST_THROW_EXCEPTION(LogicException(oss.str()));
    }

    const pos_t centerPos( isBp1 ? centerPos1 : centerPos2 );

    const pos_t insertSize(sv.insertSeq.size());

    // total impact of the alt allele on template size:
    const pos_t altShift((centerPos2-centerPos1)-insertSize);

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
        const pos_t minFrag(static_cast<pos_t>(pRange.min));
        const pos_t maxFrag(static_cast<pos_t>(pRange.max));

        const pos_t maxSupportedFrag(maxFrag-pairOpt.minFragSupport);

        const pos_t beginPos(centerPos-maxSupportedFrag);
        const pos_t endPos(centerPos+maxSupportedFrag+1);
#ifdef DEBUG_MEGAPAIR
        log_os << __FUNCTION__ << ": pair scan begin/end: " << beginPos << " " << endPos << "\n";
#endif

        /// This could occur if the fragment distribution is incredibly small --
        /// we effectively can't make use of pairs in this case:
        if (beginPos >= endPos) continue;

        const SizeDistribution& fragDistro(_readScanner.getFragSizeDistro(bamIndex));

        // set bam stream to new search interval:
        bamStream.set_new_region(sv.bp1.interval.tid, beginPos, endPos);

        while (bamStream.next())
        {
            const bam_record& bamRead(*(bamStream.get_record_ptr()));

            if (bamRead.is_filter()) continue;
            if (bamRead.is_dup()) continue;
            if (bamRead.is_secondary()) continue;
            if (bamRead.is_supplement()) continue;

            if (bamRead.is_unmapped() || bamRead.is_mate_unmapped()) continue;

            /// check for standard innie orientation:
            if (! is_innie_pair(bamRead)) continue;

#ifdef DEBUG_MEGAPAIR
            log_os << __FUNCTION__ << ": read: " << bamRead << "\n";
#endif

            /// check if fragment is too big or too small:
            const int templateSize(std::abs(bamRead.template_size()));
            const int altTemplateSize(templateSize-altShift);
            if (altTemplateSize < minFrag) continue;
            if (altTemplateSize > maxFrag) continue;

            // count only from the down stream reads
            const bool isFirstBamRead(isFirstRead(bamRead));

            // get fragment range:
            pos_t fragBeginRefPos(0);
            if (isFirstBamRead)
            {
                fragBeginRefPos=bamRead.pos()-1;
            }
            else
            {
                fragBeginRefPos=bamRead.mate_pos()-1;
            }

            const pos_t fragEndRefPos(fragBeginRefPos+templateSize);

            if (fragBeginRefPos > fragEndRefPos)
            {
                std::ostringstream oss;
                oss << "ERROR: Failed to parse fragment range from bam record. Frag begin,end: " << fragBeginRefPos << " " << fragEndRefPos << " bamRecord: " << bamRead << "\n";
                BOOST_THROW_EXCEPTION(LogicException(oss.str()));
            }

            {
                const pos_t fragOverlap(std::min((1+centerPos1-fragBeginRefPos), (fragEndRefPos-centerPos2)));
#ifdef DEBUG_MEGAPAIR
                log_os << __FUNCTION__ << ": frag begin/end/overlap: " << fragBeginRefPos << " " << fragEndRefPos << " " << fragOverlap << "\n";
#endif
                if (fragOverlap < pairOpt.minFragSupport) continue;
            }

            SVFragmentEvidence& fragment(evidence.getSample(isTumor)[bamRead.qname()]);

            SVFragmentEvidenceRead& evRead(fragment.getRead(bamRead.is_first()));
            setReadEvidence(minMapQ, bamRead, evRead);

            setAlleleFrag(fragDistro, altTemplateSize, fragment.alt.getBp(isBp1));

            // when an alt entry is made for a fragment, we /*always*/ create correponding ref entry
            // in theory this will get picked up by the ref scanner anyway, but the cost of missing this
            // is all sorts of really bad somatic FNs
            setAlleleFrag(fragDistro, templateSize, fragment.ref.getBp(isBp1));

            if (! isFirstBamRead) continue;
            if (! evRead.isAnchored) continue;

            /// old tracker:
            if (isBp1)
            {
                sample.alt.bp1SpanReadCount++;
            }
            else
            {
                sample.alt.bp2SpanReadCount++;
            }
        }
    }
}



/// get reference allele pair support at a single breakend:
///
void
SVScorer::
getSVRefPairSupport(
    const PairOptions& pairOpt,
    const SVBreakend& bp,
    const bool isBp1,
    SVScoreInfo& baseInfo,
    SVEvidence& evidence)
{
    /// search for all read pairs supporting the reference allele
    ///
    /// APPROXIMATION: for imprecise and precise variants treat the breakend locations as the center of the
    ///  breakend interval.
    ///
    /// TODO: improve on the approx above
    ///
    const pos_t centerPos(bp.interval.range.center_pos());


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
        const pos_t minFrag(pRange.min);
        const pos_t maxFrag(pRange.max);

        const SizeDistribution& fragDistro(_readScanner.getFragSizeDistro(bamIndex));

        const pos_t maxSupportedFrag(maxFrag-pairOpt.minFragSupport);

        const pos_t beginPos(centerPos-maxSupportedFrag);
        const pos_t endPos(centerPos+maxSupportedFrag+1);
#ifdef DEBUG_MEGAPAIR
        log_os << __FUNCTION__ << ": pair scan begin/end: " << beginPos << " " << endPos << "\n";
#endif

        /// This could occur if the fragment distribution is incredibly small --
        /// we effectively can't make use of pairs in this case:
        if (beginPos >= endPos) continue;

        // set bam stream to new search interval:
        bamStream.set_new_region(bp.interval.tid, beginPos, endPos);

        while (bamStream.next())
        {
            const bam_record& bamRead(*(bamStream.get_record_ptr()));

            if (bamRead.is_filter()) continue;
            if (bamRead.is_dup()) continue;
            if (bamRead.is_secondary()) continue;
            if (bamRead.is_supplement()) continue;

            if (bamRead.is_unmapped() || bamRead.is_mate_unmapped()) continue;

            /// check for standard innie orientation:
            if (! is_innie_pair(bamRead)) continue;

#ifdef DEBUG_MEGAPAIR
            log_os << __FUNCTION__ << ": read: " << bamRead << "\n";
#endif

            /// check if fragment is too big or too small:
            const int templateSize(std::abs(bamRead.template_size()));
            if (templateSize < minFrag) continue;
            if (templateSize > maxFrag) continue;

            // count only from the down stream read unless the mate-pos goes past center-pos
            const bool isLeftMost(bamRead.pos() < bamRead.mate_pos());
            const bool isRead1Tie((bamRead.pos() == bamRead.mate_pos()) && bamRead.is_first());
            const bool isDefaultSelected(isLeftMost || isRead1Tie);

            const bool isMateBeforeCenter(bamRead.mate_pos() < centerPos);

            bool isDoubleCountSkip(false);
            if ( isDefaultSelected && isMateBeforeCenter ) isDoubleCountSkip=true;
            if ( (!isDefaultSelected) && (!isMateBeforeCenter) ) isDoubleCountSkip=true;

            // get fragment range:
            pos_t fragBeginRefPos(0);
            if (isLeftMost)
            {
                fragBeginRefPos=bamRead.pos()-1;
            }
            else
            {
                fragBeginRefPos=bamRead.mate_pos()-1;
            }

            const pos_t fragEndRefPos(fragBeginRefPos+templateSize);

            if (fragBeginRefPos > fragEndRefPos)
            {
                using namespace illumina::common;

                std::ostringstream oss;
                oss << "ERROR: Failed to parse fragment range from bam record. Frag begin,end: " << fragBeginRefPos << " " << fragEndRefPos << " bamRecord: " << bamRead << "\n";
                BOOST_THROW_EXCEPTION(LogicException(oss.str()));
            }

            {
                const pos_t fragOverlap(std::min((1+centerPos-fragBeginRefPos), (fragEndRefPos-centerPos)));
#ifdef DEBUG_MEGAPAIR
                log_os << __FUNCTION__ << ": frag begin/end/overlap: " << fragBeginRefPos << " " << fragEndRefPos << " " << fragOverlap << "\n";
#endif
                if (fragOverlap < pairOpt.minFragSupport) continue;
            }

            SVFragmentEvidence& fragment(evidence.getSample(isTumor)[bamRead.qname()]);

            SVFragmentEvidenceRead& evRead(fragment.getRead(bamRead.is_first()));
            setReadEvidence(minMapQ, bamRead, evRead);

            setAlleleFrag(fragDistro, templateSize, fragment.ref.getBp(isBp1));

            if (isDoubleCountSkip) continue;
            if (! evRead.isAnchored) continue;

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
finishAllelePairSupport(
    SVSampleAlleleInfo& allele)
{
    allele.spanPairCount = std::min(allele.bp1SpanReadCount, allele.bp2SpanReadCount);
}



void
SVScorer::
getSVAltPairSupport(
    const PairOptions& pairOpt,
    const SVCandidate& sv,
    SVScoreInfo& baseInfo,
    SVEvidence& evidence)
{
    getSimpleSVAltPairSupport(pairOpt, sv, true, baseInfo, evidence);
    getSimpleSVAltPairSupport(pairOpt, sv, false, baseInfo, evidence);

    finishAllelePairSupport(baseInfo.tumor.alt);
    finishAllelePairSupport(baseInfo.normal.alt);
}



void
SVScorer::
getSVRefPairSupport(
    const PairOptions& pairOpt,
    const SVCandidate& sv,
    SVScoreInfo& baseInfo,
    SVEvidence& evidence)
{
    getSVRefPairSupport(pairOpt, sv.bp1, true, baseInfo, evidence);
    getSVRefPairSupport(pairOpt, sv.bp2, false, baseInfo, evidence);

    finishAllelePairSupport(baseInfo.tumor.ref);
    finishAllelePairSupport(baseInfo.normal.ref);
}


struct SpanReadInfo
{
    SpanReadInfo() :
        isFwdStrand(true),
        readSize(0)
    {}

    GenomeInterval interval;
    bool isFwdStrand;
    unsigned readSize;
};



static
void
getFragInfo(
    const bam_record& localRead,
    SpanReadInfo& local,
    SpanReadInfo& remote)
{
    using namespace ALIGNPATH;

    // local read:
    local.isFwdStrand = localRead.is_fwd_strand();
    local.readSize = localRead.read_size();
    local.interval.tid = localRead.target_id();
    const pos_t localBeginPos(localRead.pos()-1);

    // get cigar:
    path_t localPath;
    bam_cigar_to_apath(localRead.raw_cigar(), localRead.n_cigar(), localPath);

    const pos_t localEndPos(localBeginPos+apath_ref_length(localPath));

    local.interval.range.set_range(localBeginPos,localEndPos);

    // remote read:
    remote.isFwdStrand = localRead.is_mate_fwd_strand();
    remote.readSize = local.readSize;
    remote.interval.tid = localRead.mate_target_id();
    const pos_t remoteBeginPos(localRead.mate_pos()-1);

    // approximate end-point of remote read:
    const pos_t remoteEndPos(remoteBeginPos+localRead.read_size());

    remote.interval.range.set_range(remoteBeginPos,remoteEndPos);
}



static
void
getFragInfo(
    const SVCandidateSetReadPair& pair,
    SpanReadInfo& read1,
    SpanReadInfo& read2)
{
    using namespace ALIGNPATH;

    if (pair.read1.isSet())
    {
        getFragInfo(pair.read1.bamrec, read1, read2);

        if (pair.read2.isSet())
        {
            const bam_record& bamRead2(pair.read2.bamrec);

            read2.readSize = bamRead2.read_size();

            // get cigar:
            path_t apath2;
            bam_cigar_to_apath(bamRead2.raw_cigar(), bamRead2.n_cigar(), apath2);

            read2.interval.range.set_end_pos(read2.interval.range.begin_pos() + apath_ref_length(apath2));
        }
    }
    else if (pair.read2.isSet())
    {
        getFragInfo(pair.read2.bamrec, read2, read1);
    }
    else
    {
        assert(false && "Neither fragment read found");
    }
}



struct SpanTerminal
{
    SpanTerminal() :
        tid(0),
        pos(0),
        isFwdStrand(true),
        readSize(0)
    {}

    int32_t tid;
    pos_t pos;
    bool isFwdStrand;
    unsigned readSize;
};


#ifdef DEBUG_PAIR
static
std::ostream&
operator<<(std::ostream& os, const SpanTerminal& st)
{
    os << "tid: " << st.tid
       << " pos: "<< st.pos
       << " isFwdStrand: " << st.isFwdStrand
       << " readSize: " << st.readSize;
    return os;
}
#endif



static
void
getTerminal(
    const SpanReadInfo& rinfo,
    SpanTerminal& fterm)
{
    fterm.tid = rinfo.interval.tid;
    fterm.isFwdStrand = rinfo.isFwdStrand;
    fterm.pos = ( fterm.isFwdStrand ? rinfo.interval.range.begin_pos() : rinfo.interval.range.end_pos() );
    fterm.readSize = rinfo.readSize;
}



/// double check that a read-pair supports an sv, and if so what is the fragment length prob?
static
void
getFragProb(
    const PairOptions& pairOpt,
    const SVCandidate& sv,
    const SVCandidateSetReadPair& pair,
    const SizeDistribution& fragDistro,
    const bool isStrictMatch,
    bool& isFragSupportSV,
    float& fragProb)
{
#ifdef DEBUG_PAIR
    static const std::string logtag("getFragProb: ");
#endif

    isFragSupportSV=false;
    fragProb=0.;

    SpanReadInfo read1;
    SpanReadInfo read2;
    getFragInfo(pair, read1, read2);

    // define the end-points of fragment:
    SpanTerminal frag1;
    getTerminal(read1,frag1);

    SpanTerminal frag2;
    getTerminal(read2,frag2);

    const pos_t bp1pos(sv.bp1.interval.range.center_pos());
    const pos_t bp2pos(sv.bp2.interval.range.center_pos());

    // match bp to frag
    bool isBpFragReversed(false);

    if (frag1.tid != sv.bp1.interval.tid)
    {
        isBpFragReversed=true;
    }
    else if (frag1.isFwdStrand != (sv.bp1.state == SVBreakendState::RIGHT_OPEN) )
    {
        isBpFragReversed=true;
    }
    else if (frag1.isFwdStrand == frag2.isFwdStrand)
    {
        if ((frag1.pos < frag2.pos) != (bp1pos < bp2pos))
        {
            if (frag1.pos != frag2.pos)
            {
                isBpFragReversed=true;
            }
        }
    }

    if (isBpFragReversed)
    {
        std::swap(frag1,frag2);
#ifdef DEBUG_PAIR
        log_os << logtag << "swapping fragments\n";
#endif
    }

#ifdef DEBUG_PAIR
    log_os << logtag << "pair: " << pair << "\n";
    log_os << logtag << "sv: " << sv << "\n";
    log_os << logtag << "frag1: " << frag1 << "\n";
    log_os << logtag << "frag2: " << frag2 << "\n";
#endif

    // QC the frag/bp matchup:
    {
        std::string errorMsg;
        if (frag1.tid != frag2.tid)
        {
            if (frag1.tid != sv.bp1.interval.tid)
            {
                errorMsg = "Can't match evidence read chrom to sv-candidate bp1.";
            }
            if (frag2.tid != sv.bp2.interval.tid)
            {
                errorMsg = "Can't match evidence read chrom to sv-candidate bp2.";
            }
        }
        else if (frag1.isFwdStrand != frag2.isFwdStrand)
        {
            if ( frag1.isFwdStrand != (sv.bp1.state == SVBreakendState::RIGHT_OPEN) )
            {
                errorMsg = "Can't match evidence read strand to sv-candidate bp1";
            }
            if ( frag2.isFwdStrand != (sv.bp2.state == SVBreakendState::RIGHT_OPEN) )
            {
                errorMsg = "Can't match evidence read strand to sv-candidate bp2";
            }
        }
        else
        {
            if ( (frag1.pos < frag2.pos) != (bp1pos < bp2pos) )
            {
                if (frag1.pos != frag2.pos)
                {
                    errorMsg = "Can't match read pair positions to sv-candidate.";
                }
            }
        }

        if (! errorMsg.empty())
        {
            if (! isStrictMatch) return;

            using namespace illumina::common;

            std::ostringstream oss;
            oss << "ERROR: " << errorMsg  << "\n"
                << "\tcandidate-sv: " << sv
                << "\tread-pair: " << pair
                << "\n";
            BOOST_THROW_EXCEPTION(LogicException(oss.str()));
        }
    }

    pos_t frag1Size(bp1pos-frag1.pos);
    if (! frag1.isFwdStrand) frag1Size *= -1;

    pos_t frag2Size(bp2pos-frag2.pos);
    if (! frag2.isFwdStrand) frag2Size *= -1;

#ifdef DEBUG_PAIR
    log_os << logtag << "frag1size,frag2size: " << frag1Size << " " << frag2Size << "\n";
#endif

    if (frag1Size < std::min(pairOpt.minFragSupport, static_cast<pos_t>(frag1.readSize))) return;
    if (frag2Size < std::min(pairOpt.minFragSupport, static_cast<pos_t>(frag2.readSize))) return;


    fragProb=fragDistro.cdf(frag1Size+frag2Size);
#ifdef DEBUG_PAIR
    log_os << logtag << "cdf: " << fragProb << " final: " << std::min(fragProb, (1-fragProb)) << "\n";
#endif
    fragProb = std::min(fragProb, (1-fragProb));

    /// TODO: any cases where fragProb is 0 should be some soft of mulit-SV error artifact (like a large CIGAR indel in one of the reads of the pair)
    ///     try to improve this case -- ideally we can account for such events.
    if (fragProb > 0.)
    {
        isFragSupportSV = true;
    }
}



// count the read pairs supporting the alternate allele in each sample, using data we already produced during candidate generation:
//
void
SVScorer::
processExistingAltPairInfo(
    const PairOptions& pairOpt,
    const SVCandidateSetData& svData,
    const SVCandidate& sv,
    SVScoreInfo& baseInfo,
    SVEvidence& evidence)
{
    const unsigned minMapQ(_readScanner.getMinMapQ());

    const unsigned bamCount(_bamStreams.size());
    for (unsigned bamIndex(0); bamIndex < bamCount; ++bamIndex)
    {
        const bool isTumor(_isAlignmentTumor[bamIndex]);
        SVSampleInfo& sample(isTumor ? baseInfo.tumor : baseInfo.normal);

        const SizeDistribution& fragDistro(_readScanner.getFragSizeDistro(bamIndex));

        const SVCandidateSetReadPairSampleGroup& svDataGroup(svData.getDataGroup(bamIndex));
        BOOST_FOREACH(const SVCandidateSetReadPair& pair, svDataGroup)
        {
            // is this read pair associated with this candidateIndex? (each read pair can be associated with multiple candidates)
            unsigned linkIndex(0);
            {
                bool isIndexFound(false);
                BOOST_FOREACH(const SVPairAssociation& sva, pair.svLink)
                {
                    if (sv.candidateIndex == sva.index)
                    {
                        isIndexFound = true;
                        break;
                    }
                    linkIndex++;
                }

                if (! isIndexFound) continue;
            }

            /// do we assert (strict) or skip (non-strict) on non-matching pair/sv associations?
            bool isStrictMatch(true);
            if (pair.svLink.size() > 1)
            {
                if (! SVEvidenceType::isPairType(pair.svLink[linkIndex].evtype)) isStrictMatch=false;
            }

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

#ifdef DEBUG_PAIR
            static const std::string logtag("processExistingAltPairInfo: ");
            log_os << logtag << "Finding alt pair evidence for svIndex: " << sv.candidateIndex << "  qname: " << qname << "\n";
#endif

            SVFragmentEvidence& fragment(evidence.getSample(isTumor)[qname]);
            SVFragmentEvidenceAllele& alt(fragment.alt);

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

            /// get fragment prob, and possibly withdraw fragment support based on refined sv breakend coordinates:
            bool isFragSupportSV(false);
            float fragProb(0);
            if (pair.read1.isSet() && pair.read2.isSet())
            {
                getFragProb(pairOpt, sv, pair, fragDistro, isStrictMatch, isFragSupportSV, fragProb);
            }

            if (! isFragSupportSV) continue;

            /// TODO: if fragProb is zero this should be a bug -- follow-up to see if we can make this an assert(fragProb > 0.) instead
            if (fragProb <= 0.) continue;

            // for all large spanning events -- we don't test for pair support of the two breakends separately -- this could be
            // beneficial if there was an unusually large insertion associated with the event. For now we approximate that
            // these events will mostly not have very large insertions.
            //
            alt.bp1.isFragmentSupport = true;
            alt.bp1.fragLengthProb = fragProb;

            alt.bp2.isFragmentSupport = true;
            alt.bp2.fragLengthProb = fragProb;
        }
    }
}



void
SVScorer::
getSVPairSupport(
    const SVCandidateSetData& svData,
    const SVCandidateAssemblyData& assemblyData,
    const SVCandidate& sv,
    SVScoreInfo& baseInfo,
    SVEvidence& evidence)
{
    static const PairOptions pairOpt;

#ifdef DEBUG_PAIR
    static const std::string logtag("getSVPairSupport: ");
    log_os << logtag << "starting alt pair search for sv: " << sv << "\n";
#endif

    if (svData.isSkipped()) return;

    if (assemblyData.isCandidateSpanning)
    {
        // count the read pairs supporting the alternate allele in each sample
        // using data we already produced during candidate generation:
        //
        processExistingAltPairInfo(pairOpt, svData, sv, baseInfo, evidence);
    }
    else
    {
        // for SVs which were assembled without a pair-driven prior hypothesis,
        // we need to go back to the bam and and find any supporting alt read-pairs
        getSVAltPairSupport(pairOpt, sv, baseInfo, evidence);
    }

    // count the read pairs supporting the reference allele on each breakend in each sample:
    //
    getSVRefPairSupport(pairOpt, sv, baseInfo, evidence);
}

