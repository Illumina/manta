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
#include "SVScorerShared.hh"
#include "SVScorePairAltProcessor.hh"
#include "SVScorePairRefProcessor.hh"

#include "blt_util/align_path_bam_util.hh"
#include "blt_util/bam_streamer.hh"
#include "blt_util/bam_record_util.hh"
#include "common/Exceptions.hh"
#include "manta/SVCandidateUtil.hh"
#include "svgraph/GenomeIntervalUtil.hh"

#include "boost/foreach.hpp"
#include "boost/make_shared.hpp"

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
processBamProcList(
    const std::vector<SVScorer::streamPtr>& bamList,
    std::vector<SVScorer::pairProcPtr>& pairProcList)
{
    const unsigned bamCount(bamList.size());
    const unsigned bamProcCount(pairProcList.size());

    for (unsigned bamIndex(0); bamIndex < bamCount; ++bamIndex)
    {
        // get the minimum set of scan intervals (this should almost always be 1!)
        std::vector<GenomeInterval> scanIntervals;
        std::vector<unsigned> intervalMap;
        {
            BOOST_FOREACH(SVScorer::pairProcPtr& bpp, pairProcList)
            {
                const GenomeInterval& interval(bpp->nextBamIndex(bamIndex));
                if (interval.range.size() < 1) continue;

                scanIntervals.push_back(interval);
            }

            intervalMap = intervalCompressor(scanIntervals);
        }

        bam_streamer& bamStream(*bamList[bamIndex]);

        const unsigned intervalCount(scanIntervals.size());
        for (unsigned intervalIndex(0); intervalIndex<intervalCount; ++intervalIndex)
        {
            const GenomeInterval& scanInterval(scanIntervals[intervalIndex]);
            if (scanInterval.range.begin_pos() >= scanInterval.range.end_pos()) continue;

            // set bam stream to new search interval:
            bamStream.set_new_region(scanInterval.tid, scanInterval.range.begin_pos(), scanInterval.range.end_pos());

            /// define the procs where' going to handle in this interval:
            std::vector<unsigned> targetProcs;
            for (unsigned procIndex(0); procIndex<bamProcCount; ++procIndex)
            {
                if (intervalMap[procIndex] == intervalIndex)
                {
                    targetProcs.push_back(procIndex);
                }
            }

            while (bamStream.next())
            {
                const bam_record& bamRead(*(bamStream.get_record_ptr()));

                /// this filter is common to all targetProcs:
                if (SVScorer::pairProcPtr::value_type::isSkipRecordCore(bamRead)) continue;

                BOOST_FOREACH(const unsigned procIndex, targetProcs)
                {
                    SVScorer::pairProcPtr& bpp(pairProcList[procIndex]);

                    if (bpp->isSkipRecord(bamRead)) continue;
                    bpp->processClearedRecord(bamRead);
                }
            }
        }
    }
}



void
SVScorer::
getSVAltPairSupport(
    const PairOptions& pairOpt,
    const SVCandidateAssemblyData& assemblyData,
    const SVCandidate& sv,
    SVEvidence& evidence,
    std::vector<pairProcPtr>& pairProcList)
{
    pairProcPtr bp1Ptr(new SVScorePairAltProcessor(_scanOpt, _refineOpt, _isAlignmentTumor, _readScanner, pairOpt, assemblyData, sv, true, evidence));
    pairProcPtr bp2Ptr(new SVScorePairAltProcessor(_scanOpt, _refineOpt, _isAlignmentTumor, _readScanner, pairOpt, assemblyData, sv, false, evidence));

    pairProcList.push_back(bp1Ptr);
    pairProcList.push_back(bp2Ptr);
}



void
SVScorer::
getSVRefPairSupport(
    const PairOptions& pairOpt,
    const SVCandidate& sv,
    SVEvidence& evidence,
    std::vector<pairProcPtr>& pairProcList)
{
    pairProcPtr bp1Ptr(new SVScorePairRefProcessor(_isAlignmentTumor, _readScanner, pairOpt, sv, true, evidence));
    pairProcPtr bp2Ptr(new SVScorePairRefProcessor(_isAlignmentTumor, _readScanner, pairOpt, sv, false, evidence));

    pairProcList.push_back(bp1Ptr);
    pairProcList.push_back(bp2Ptr);
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



/// fill in SpanReadInfo as accurately as possible depending on
/// whether one or both of the read pair's bam records have been found:
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



/// read pairs are abstracted to two terminals for the purpose of
/// fragment size estimation in the context of the alternate allele:
/// tid+pos represent one of the two extreme ends of the fragment in
/// genomic chromosome+position coordinates
///
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



/// convert SpanReadInfo to SpanTerminal
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
            oss << "ERROR: " << errorMsg  << '\n'
                << "\tcandidate-sv: " << sv
                << "\tread-pair: " << pair
                << '\n';
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

    if (frag1Size < pairOpt.minFragSupport) return;
    if (frag2Size < pairOpt.minFragSupport) return;

    fragProb=fragDistro.cdf(frag1Size+frag2Size);
#ifdef DEBUG_PAIR
    log_os << logtag << "cdf: " << fragProb << " final: " << std::min(fragProb, (1-fragProb)) << "\n";
#endif
    fragProb = std::min(fragProb, (1-fragProb));

    /// TODO: any cases where fragProb is 0 or extremely small should be some
    /// sort of mulit-SV error artifact (like a large CIGAR indel in one of the
    /// reads of the pair) try to improve this case -- ideally we can account
    /// for such events.
    if (fragProb >= pairOpt.minFragProb)
    {
        isFragSupportSV = true;
    }

#ifdef DEBUG_PAIR
    log_os << logtag << "isSupportSV: " << isFragSupportSV << "\n";
#endif
}



// count the read pairs supporting the alternate allele in each sample, using data we already produced during candidate generation:
//
void
SVScorer::
processExistingAltPairInfo(
    const PairOptions& pairOpt,
    const SVCandidateSetData& svData,
    const SVCandidate& sv,
    SVEvidence& evidence)
{
    const unsigned minMapQ(_readScanner.getMinMapQ());
    const unsigned minTier2MapQ(_readScanner.getMinTier2MapQ());

    const unsigned bamCount(_bamStreams.size());
    for (unsigned bamIndex(0); bamIndex < bamCount; ++bamIndex)
    {
        const bool isTumor(_isAlignmentTumor[bamIndex]);

        const SizeDistribution& fragDistro(_readScanner.getFragSizeDistro(bamIndex));

        const SVCandidateSetReadPairSampleGroup& svDataGroup(svData.getDataGroup(bamIndex));
        BOOST_FOREACH(const SVCandidateSetReadPair& pair, svDataGroup)
        {
            /// at least one read of the pair must have been found:
            assert(pair.read1.isSet() || pair.read2.isSet());

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
            assert(pair.svLink.size() > linkIndex);

            const bool isPairType(SVEvidenceType::isPairType(pair.svLink[linkIndex].evtype));

            /// if the evidence comes from a read pair observation, a very strict matching criteria
            /// is enforced between this pair and the SV candidate. If the read pair association comes from
            /// a CIGAR string for instance, the pair will not necessarily support the candidate
            ///
            const bool isStrictMatch(isPairType);

            const std::string qname(pair.qname());

#ifdef DEBUG_PAIR
            static const std::string logtag("processExistingAltPairInfo: ");
            log_os << logtag << "Finding alt pair evidence for svIndex: " << sv.candidateIndex << "  qname: " << qname << "\n";
#endif

            SVFragmentEvidence& fragment(evidence.getSample(isTumor)[qname]);
            SVFragmentEvidenceAllele& alt(fragment.alt);

            static const bool isShadow(false);
            if (pair.read1.isSet())
            {
                setReadEvidence(minMapQ, minTier2MapQ, pair.read1.bamrec, isShadow, fragment.read1);
            }

            if (pair.read2.isSet())
            {
                setReadEvidence(minMapQ, minTier2MapQ, pair.read2.bamrec, isShadow, fragment.read2);
            }

            /// get fragment prob, and possibly withdraw fragment support based on refined sv breakend coordinates:
            bool isFragSupportSV(false);
            float fragProb(0);
            getFragProb(pairOpt, sv, pair, fragDistro, isStrictMatch, isFragSupportSV, fragProb);

            if (! isFragSupportSV) continue;

            /// TODO: if fragProb is zero this should be a bug -- follow-up to see if we can make this an assert(fragProb > 0.) instead
            if (fragProb <= 0.)
            {
#ifdef DEBUG_PAIR
                log_os << logtag << "Fragment with fragProb=0! " << sv.candidateIndex << "  qname: " << qname << "\n";
#endif
                continue;
            }

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
    SVEvidence& evidence)
{
    const PairOptions pairOpt(_isRNA);

#ifdef DEBUG_PAIR
    static const std::string logtag("getSVPairSupport: ");
    log_os << logtag << "starting alt pair search for sv: " << sv << "\n";
#endif

    std::vector<pairProcPtr> pairProcList;

    if (assemblyData.isCandidateSpanning)
    {
        // count the read pairs supporting the alternate allele in each sample
        // using data we already produced during candidate generation:
        //
        processExistingAltPairInfo(pairOpt, svData, sv, evidence);
    }
    else
    {
        // for SVs which were assembled without a pair-driven prior hypothesis,
        // we need to go back to the bam and and find any supporting alt read-pairs
        getSVAltPairSupport(pairOpt, assemblyData, sv, evidence, pairProcList);
    }

    // count the read pairs supporting the reference allele on each breakend in each sample:
    //
    getSVRefPairSupport(pairOpt, sv, evidence, pairProcList);

    // execute bam scanning for all pairs:
    //
    processBamProcList(_bamStreams, pairProcList);
}
