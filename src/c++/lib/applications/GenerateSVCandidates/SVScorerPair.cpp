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
#include "boost/foreach.hpp"

#include <iostream>
#include <sstream>
#include <string>


//#define DEBUG_SVS

#ifdef DEBUG_SVS
#include "blt_util/log.hh"
#endif


//#define DEBUG_PAIR

#ifdef DEBUG_PAIR
#include "blt_util/log.hh"
#endif



/// get reference allele support at a single breakend:
///
void
SVScorer::
getSVRefPairSupport(
    const PairOptions& pairOpt,
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

        const SizeDistribution& fragDistro(_readScanner.getFragSizeDistro(bamIndex));

        const uint32_t beginPos(centerPos-maxFrag+pairOpt.minFragSupport);
        const uint32_t endPos(centerPos+maxFrag-pairOpt.minFragSupport+1);

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

            const unsigned fragLength(std::abs(bamRead.template_size()));
            const pos_t fragEnd(fragBegin+fragLength);

            if (fragBegin > fragEnd)
            {
                using namespace illumina::common;

                std::ostringstream oss;
                oss << "ERROR: Failed to parse fragment range from bam record. Frag begin,end: " << fragBegin << " " << fragEnd << " bamRecord: " << bamRead << "\n";
                BOOST_THROW_EXCEPTION(LogicException(oss.str()));
            }

            {
                const pos_t fragOverlap(std::min((centerPos-fragBegin), (fragEnd-centerPos)));
                if (fragOverlap < std::min(static_cast<pos_t>(bamRead.read_size()), pairOpt.minFragSupport)) continue;
            }

            SVFragmentEvidence& fragment(evidence.getSample(isTumor)[bamRead.qname()]);
            SVFragmentEvidenceAllele& ref(fragment.ref);

            setReadEvidence(minMapQ, bamRead, fragment.getRead(bamRead.is_first()));

            {
                float fragProb(fragDistro.cdf(fragLength));
                fragProb = std::min(fragProb, (1-fragProb));

                if (isBp1)
                {
                    ref.bp1.isFragmentSupport = true;
                    ref.bp1.fragLengthProb = fragProb;
                }
                else
                {
                    ref.bp2.isFragmentSupport = true;
                    ref.bp2.fragLengthProb = fragProb;
                }
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
    const PairOptions& pairOpt,
    const SVCandidate& sv,
    SVScoreInfo& baseInfo,
    SVEvidence& evidence)
{
    getSVRefPairSupport(pairOpt, sv.bp1, baseInfo, evidence, true);
    getSVRefPairSupport(pairOpt, sv.bp2, baseInfo, evidence, false);

    finishSamplePairSupport(baseInfo.tumor);
    finishSamplePairSupport(baseInfo.normal);
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
    else if(pair.read2.isSet())
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
        isFwd(true),
        readSize(0)
    {}

    int32_t tid;
    pos_t pos;
    bool isFwd;
    unsigned readSize;
};


std::ostream&
operator<<(std::ostream& os, const SpanTerminal& st)
{
    os << "tid: " << st.tid << " pos: " << st.pos << " isF: " << st.isFwd << " readSize: " << st.readSize;
    return os;
}



static
void
getTerminal(
    const SpanReadInfo& rinfo,
    SpanTerminal& fterm)
{
    fterm.tid = rinfo.interval.tid;
    fterm.isFwd = rinfo.isFwdStrand;
    fterm.pos = ( fterm.isFwd ? rinfo.interval.range.begin_pos() : rinfo.interval.range.end_pos() );
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
    bool& isFragSupportSV,
    float& fragProb)
{
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
    else if(frag1.isFwd != (sv.bp1.state == SVBreakendState::RIGHT_OPEN) )
    {
        isBpFragReversed=true;
    }
    else if((frag1.pos < frag2.pos) != (bp1pos < bp2pos))
    {
        isBpFragReversed=true;
    }

    if (isBpFragReversed)
    {
        std::swap(frag1,frag2);
    }

#ifdef DEBUG_PAIR
    static const std::string logtag("getFragProb: ");
    log_os << logtag << "read1: ";
    if (pair.read1.isSet())
    {
        log_os << pair.read1.bamrec << "\n";
    }
    else
    {
        log_os << "UNKNOWN\n";
    }

    log_os << logtag << "read2: ";
    if (pair.read1.isSet())
    {
        log_os << pair.read2.bamrec << "\n";
    }
    else
    {
        log_os << "UNKNOWN\n";
    }
    log_os << logtag << "sv: " << sv << "\n";
    log_os << logtag << "frag1: " << frag1 << "\n";
    log_os << logtag << "frag2: " << frag2 << "\n";
#endif

    // QC the frag/bp matchup:
    if (frag1.tid != frag2.tid)
    {
        /// TODO:: we should be able to assert this condition with no return.. there's an occational bad read/sv matchup
        if (frag1.tid != sv.bp1.interval.tid) return;
        if (frag2.tid != sv.bp2.interval.tid) return;
    }
    else if(frag1.isFwd != frag2.isFwd)
    {
        /// TODO:: we should be able to assert this condition with no return.. there's an occational bad read/sv matchup
        if ( frag1.isFwd != (sv.bp1.state == SVBreakendState::RIGHT_OPEN) ) return;
        if ( frag2.isFwd != (sv.bp2.state == SVBreakendState::RIGHT_OPEN) ) return;;
    }
    else
    {
        assert ( (frag1.pos < frag2.pos) == (bp1pos < bp2pos) );
    }


    pos_t frag1Size(bp1pos-frag1.pos);
    if (! frag1.isFwd) frag1Size *= -1;

    pos_t frag2Size(bp2pos-frag2.pos);
    if (! frag2.isFwd) frag2Size *= -1;

#ifdef DEBUG_PAIR
    log_os << logtag << "frag1size,frag2size: " << frag1Size << " " << frag2Size << "\n";
#endif

    if (frag1Size < std::min(pairOpt.minFragSupport, static_cast<pos_t>(frag1.readSize))) return;
    if (frag2Size < std::min(pairOpt.minFragSupport, static_cast<pos_t>(frag2.readSize))) return;


    isFragSupportSV = true;

    fragProb=fragDistro.cdf(frag1Size+frag2Size);
#ifdef DEBUG_PAIR
    log_os << logtag << "cdf: " << fragProb << " final: " << std::min(fragProb, (1-fragProb)) << "\n";
#endif
    fragProb = std::min(fragProb, (1-fragProb));
}




void
SVScorer::
getSVPairSupport(
    const SVCandidateSetData& svData,
    const SVCandidate& sv,
    SVScoreInfo& baseInfo,
    SVEvidence& evidence)
{
    static PairOptions pairOpt;

    const unsigned minMapQ(_readScanner.getMinMapQ());

    // count the read pairs supporting the alternate allele in each sample, using data we already produced during candidate generation:
    //
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

#ifdef DEBUG_PAIR
            log_os << "getSVPairSupport: Finding alt pair evidence for: " << qname << "\n";
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
            getFragProb(pairOpt, sv, pair, fragDistro, isFragSupportSV, fragProb);

            if (! isFragSupportSV) continue;

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

    // count the read pairs supporting the reference allele on each breakend in each sample:
    //
    getSVRefPairSupport(pairOpt, sv, baseInfo, evidence);
}
