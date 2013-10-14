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
/// \author Chris Saunders
/// \author Ole Schulz-Trieglaff
/// \author Bret Barnes
///

#include "alignment/ReadScorer.hh"
#include "blt_util/align_path_bam_util.hh"
#include "blt_util/align_path_util.hh"
#include "blt_util/bam_record_util.hh"
#include "blt_util/parse_util.hh"
#include "blt_util/string_util.hh"
#include "common/Exceptions.hh"
#include "manta/SVLocusScanner.hh"

#include "boost/foreach.hpp"


//#define DEBUG_SCANNER

//#define DEBUG_SEMI_ALIGNED
//#define DEBUG_IS_SHADOW

#ifdef DEBUG_SCANNER
#include "blt_util/log.hh"
#include <iostream>
#endif


const float SVObservationWeights::closePairFactor(4);



struct SimpleAlignment
{
    SimpleAlignment() :
        is_fwd_strand(true),
        pos(0)
    {}

    SimpleAlignment(const bam_record& bamRead) :
        is_fwd_strand(bamRead.is_fwd_strand()),
        pos(bamRead.pos()-1)
    {
        bam_cigar_to_apath(bamRead.raw_cigar(),bamRead.n_cigar(),path);
    }

    bool is_fwd_strand;
    pos_t pos;
    ALIGNPATH::path_t path;
};


struct ChromAlignment : public SimpleAlignment
{
    ChromAlignment() :
        SimpleAlignment(),
        tid(0)
    {}

    ChromAlignment(const bam_record& bamRead) :
        SimpleAlignment(bamRead),
        tid(bamRead.target_id())
    {}

    int32_t tid;
};


static
SVCandidate
GetSplitSVCandidate(
    const ReadScannerDerivOptions& dopt,
    const int32_t alignTid,
    const pos_t leftPos,
    const pos_t rightPos,
    const bool isComplex = false)
{
    SVCandidate sv;
    SVBreakend& localBreakend(sv.bp1);
    SVBreakend& remoteBreakend(sv.bp2);

    localBreakend.interval.tid = alignTid;
    remoteBreakend.interval.tid = alignTid;

    localBreakend.splitCount++;

    if (! isComplex)
    {
        remoteBreakend.splitCount++;
        localBreakend.state = SVBreakendState::RIGHT_OPEN;
        remoteBreakend.state = SVBreakendState::LEFT_OPEN;
    }
    else
    {
        localBreakend.state = SVBreakendState::COMPLEX;
        remoteBreakend.state = SVBreakendState::UNKNOWN;
    }

    localBreakend.interval.range.set_begin_pos(std::max(0,leftPos-dopt.beforeBreakend));

    if (! isComplex)
    {
        localBreakend.interval.range.set_end_pos(leftPos+dopt.afterBreakend);
    }
    else
    {
        localBreakend.interval.range.set_end_pos(rightPos+dopt.afterBreakend);
    }

    remoteBreakend.interval.range.set_begin_pos(std::max(0,rightPos-dopt.beforeBreakend));
    remoteBreakend.interval.range.set_end_pos(rightPos+dopt.afterBreakend);

    return sv;
}



/// determine if reads is open ended upstream or down stream
static
bool
isCigarUpstream(
    const ALIGNPATH::path_t& align)
{
    using namespace ALIGNPATH;
    return (apath_read_lead_size(align) > apath_read_trail_size(align));
}



static
void
updateSABreakend(
    const ReadScannerDerivOptions& dopt,
    const ChromAlignment& align,
    SVBreakend& breakend)
{
    // Need to use the match descriptors to determine if
    //  we are upstream clipped or downstream clipped. We
    //  also then need to use the strand. Below is the logic
    //
    // Upstream & +    => LEFT_OPEN
    // Upstream & -    => RIGHT_OPEN
    // DownStream & +  => RIGHT_OPEN
    // DownStream & -  => LEFT_OPEN
    //

    const bool isUpstream(isCigarUpstream(align.path));

    if (isUpstream == align.is_fwd_strand)
    {
        breakend.state = SVBreakendState::LEFT_OPEN;
    }
    else
    {
        breakend.state = SVBreakendState::RIGHT_OPEN;
    }

    breakend.interval.tid = align.tid;
    breakend.interval.range.set_begin_pos(std::max(0,align.pos-dopt.beforeBreakend));
    breakend.interval.range.set_end_pos(align.pos+dopt.afterBreakend);
}



/// get SV candidates from indels in the read alignment
static
SVCandidate
GetSplitSACandidate(
    const ReadScannerDerivOptions& dopt,
    const ChromAlignment& localAlign,
    const ChromAlignment& remoteAlign)
{
    SVCandidate sv;

    sv.bp1.splitCount++;
    sv.bp2.splitCount++;

    updateSABreakend(dopt, localAlign, sv.bp1);
    updateSABreakend(dopt, remoteAlign, sv.bp2);

    return sv;
}



typedef std::map<std::string, int32_t> chromMap_t;



static
void
getSACandidatesFromRead(
    const ReadScannerDerivOptions& dopt,
    const bam_record& localRead,
    const ChromAlignment& localAlign,
    const chromMap_t& chromToIndex,
    std::vector<SVCandidate>& candidates)
{
    using namespace ALIGNPATH;

    std::vector<std::string> saVec;
    {
        static const char satag[] = {'S','A'};
        const char* saStr(localRead.get_string_tag(satag));
        if (NULL == saStr) return;

        split_string(saStr, ';', saVec);
    }

    ChromAlignment remoteAlign;

    // For now we will only handle a single split alignment
    //  In the future we will need to sort the SA tags by order on of segements on
    //   the actual template.
    //  We may also have to toss conflicting segments that map to two different areas,
    //   or at least find some way of dealing with them.
    //std::cerr << "At: '" << saStr << "', '" << localRead << "'" << std::endl;
    //std::cerr << "Size: " << saVec.size() << std::endl;

    if (saVec.size() > 1) return;

    BOOST_FOREACH(const std::string& sa, saVec)
    {
#ifdef DEBUG_SCANNER
        log_os << "SA STRING: " << sa << "\n";
#endif
        std::vector<std::string> saDat;
        split_string(sa, ',', saDat);

        assert((saDat.size() == 6) && "Unexpected number of SA tag values");

        const chromMap_t::const_iterator ci(chromToIndex.find(saDat[0]));
        assert(ci != chromToIndex.end());

        remoteAlign.tid=(ci->second); // convert chr to int32_t via new bam header map

        remoteAlign.pos = (illumina::blt_util::parse_int_str(saDat[1])-1);

        {
            const char saStrand(saDat[2][0]); // convert to char
            assert((saStrand=='-') || (saStrand=='+'));
            remoteAlign.is_fwd_strand = (saStrand == '+');
        }

        cigar_to_apath(saDat[3].c_str(), remoteAlign.path);

        /*std::cerr << "SA Values: "
                  << ", " << saChr
                  << ", " << saPos
                  << ", " << saStrand
                  << ", " << remotePath << std::endl;
        */

        // At this point we don't care about strand
        candidates.push_back(GetSplitSACandidate(dopt, localAlign, remoteAlign));
    }
}



static
void
getSVCandidatesFromReadIndels(
    const ReadScannerOptions& opt,
    const ReadScannerDerivOptions& dopt,
    const ChromAlignment& align,
    std::vector<SVCandidate>& candidates)
{
    using namespace ALIGNPATH;

    const std::pair<unsigned,unsigned> ends(get_match_edge_segments(align.path));

    unsigned pathIndex(0);
    unsigned readOffset(0);
    pos_t refHeadPos(align.pos);

    const unsigned pathSize(align.path.size());
    while (pathIndex<pathSize)
    {
        const path_segment& ps(align.path[pathIndex]);
        const bool isBeginEdge(pathIndex<ends.first);
        const bool isEndEdge(pathIndex>ends.second);
        const bool isEdgeSegment(isBeginEdge || isEndEdge);

        // in this case, swap means combined insertion/deletion
        const bool isSwapStart(is_segment_swap_start(align.path,pathIndex));

        assert(! (isEdgeSegment && isSwapStart));

        unsigned nPathSegments(1); // number of path segments consumed
        if (isEdgeSegment)
        {
            // edge inserts are allowed for intron adjacent and grouper reads, edge deletions for intron adjacent only

            if (ps.type == INSERT)
            {
                if (ps.length >= opt.minCandidateVariantSize)
                {
                    static const bool isComplex(true);
                    candidates.push_back(GetSplitSVCandidate(dopt, align.tid, refHeadPos, refHeadPos, isComplex));
                }
            }
        }
        else if (isSwapStart)
        {
            const swap_info sinfo(align.path,pathIndex);
            if ((sinfo.delete_length >= opt.minCandidateVariantSize) || (sinfo.insert_length >= opt.minCandidateVariantSize))
            {
                candidates.push_back(GetSplitSVCandidate(dopt, align.tid, refHeadPos, refHeadPos+sinfo.delete_length));
            }

            nPathSegments = sinfo.n_seg;
        }
        else if (is_segment_type_indel(align.path[pathIndex].type))
        {
            // regular indel:

            if (ps.type == DELETE)
            {
                if (ps.length >= opt.minCandidateVariantSize)
                {
                    candidates.push_back(GetSplitSVCandidate(dopt, align.tid, refHeadPos, refHeadPos+ps.length));
                }
            }
            else if (ps.type == INSERT)
            {
                if (ps.length >= opt.minCandidateVariantSize)
                {
                    candidates.push_back(GetSplitSVCandidate(dopt, align.tid, refHeadPos, refHeadPos));
                }
            }
        }

        for (unsigned i(0); i<nPathSegments; ++i)
        {
            increment_path(align.path, pathIndex, readOffset, refHeadPos);
        }
    }
}



void
getSVBreakendCandidateClip(
    const bam_record& bamRead,
    const ALIGNPATH::path_t& apath,
    unsigned& leadingClipLen,
    unsigned& trailingClipLen,
    const uint8_t minQ,
    const float minQFrac)
{
    leadingClipLen = 0;
    trailingClipLen = 0;

    const uint8_t* qual(bamRead.qual());
    const unsigned readSize(bamRead.read_size());

    const unsigned trailingClipLenTmp(apath_soft_clip_trail_size(apath));
    if (0 != trailingClipLenTmp)
    {
        // check the quality of clipped region
        unsigned minQCount(0);
        for (unsigned pos(0); pos<trailingClipLenTmp; ++pos)
        {
            if (qual[readSize-pos-1] >= minQ) minQCount++;
        }
        if ((static_cast<float>(minQCount)/trailingClipLenTmp) >= minQFrac)
        {
            trailingClipLen = trailingClipLenTmp;
        }
    }

    const unsigned leadingClipLenTmp(apath_soft_clip_lead_size(apath));
    if (0 != leadingClipLenTmp)
    {
        // check the quality of clipped region
        unsigned minQCount(0);
        for (unsigned pos(0); pos<leadingClipLenTmp; ++pos)
        {
            if (qual[pos] >= minQ) minQCount++;
        }
        if ((static_cast<float>(minQCount)/leadingClipLenTmp) >= minQFrac)
        {
            leadingClipLen = leadingClipLenTmp;
        }
    }
}



bool
isSemiAligned(const bam_record& bamRead, const double minSemiAlignedScore)
{
    ALIGNPATH::path_t apath;
    bam_cigar_to_apath(bamRead.raw_cigar(),bamRead.n_cigar(),apath);
    const double semiAlignedScore(ReadScorer::getSemiAlignedMetric(bamRead.read_size(),apath,bamRead.qual()));
#ifdef DEBUG_SEMI_ALIGNED
    static const std::string logtag("isSemiAligned");
    log_os << logtag << " semi-aligned score=" << semiAlignedScore << " read qname=" << bamRead.qname() << " apath=" << apath <<  std::endl;
#endif
    return (semiAlignedScore>minSemiAlignedScore);
}



bool
isGoodShadow(const bam_record& bamRead,
             const uint8_t lastMapq,
             const std::string& lastQname,
             const double minSingletonMapq)
{
#ifdef DEBUG_IS_SHADOW
    static const std::string logtag("isGoodShadow");
#endif
    // shadow read should be unmapped
    if (!bamRead.is_unmapped()) return false;
    // but its partner should be aligned
    if (bamRead.is_mate_unmapped()) return false;

    static const unsigned minAvgQualShadow = 25;

    if (get_avg_quality(bamRead) < minAvgQualShadow)
    {
        return false;
    }

    if (bamRead.qname() != lastQname)
    {
        // something went wrong here, shadows should have their singleton partner
        // preceding them in the BAM file.
#ifdef DEBUG_IS_SHADOW
        log_os << logtag << " ERROR: Shadow without matching singleton : " << bamRead.qname() << " vs " << lastQname << std::endl;
#endif
        return false;
    }

    if ((unsigned int)lastMapq > minSingletonMapq)
    {
#ifdef DEBUG_IS_SHADOW
        log_os << logtag << " Found shadow!" << std::endl;
        log_os << logtag << " this mapq  = " << ((unsigned int)bamRead.map_qual())
               << " this qname = " << bamRead.qname() << std::endl;
        log_os << logtag << " last mapq  = " << ((unsigned int)lastMapq)
               << " last qname = " << lastQname << std::endl;
#endif
        return true;
    }

    return false;
}



/// get SV candidates from semi-aligned reads
static
void
getSVCandidatesFromSemiAligned(
    const ReadScannerOptions& opt,
    const bam_record& bamRead,
    const SimpleAlignment& bamAlign,
    std::vector<SVCandidate>& candidates)
{
    // semi-aligned reads don't define a full hypothesis, so they're always evidence for a 'complex' ie. undefined, event
    // in a fashion analogous to clipped reads
    static const bool isComplex(true);

    const double semiAlignedScore(ReadScorer::getSemiAlignedMetric(bamRead.read_size(),bamAlign.path,bamRead.qual()));

    //std::cout << "getSVCandidatesFromSemiAligned : semi-aligned score is " << semiAlignedScore << std::endl;
    if (semiAlignedScore>opt.minSemiAlignedScoreGraph)
    {
        const pos_t pos(bamAlign.pos);
        candidates.push_back(GetSplitSVCandidate(opt,bamRead.target_id(),pos,pos,isComplex));
    }
}



/// get SV candidates from read clipping
static
void
getSVCandidatesFromReadClip(
    const ReadScannerOptions& opt,
    const bam_record& bamRead,
    const SimpleAlignment& bamAlign,
    std::vector<SVCandidate>& candidates)
{
    unsigned leadingClipLen(0), trailingClipLen(0);
    getSVBreakendCandidateClip(bamRead, bamAlign.path, leadingClipLen, trailingClipLen);

    // soft-clipped reads don't define a full hypothesis, so they're always evidence for a 'complex' ie. undefined, event:
    static const bool isComplex(true);

    if (leadingClipLen >= opt.minSoftClipLen)
    {
        const pos_t clipPos(bamAlign.pos);
        candidates.push_back(GetSplitSVCandidate(opt,bamRead.target_id(),clipPos,clipPos,isComplex));
    }

    if (trailingClipLen >= opt.minSoftClipLen)
    {
        const pos_t clipPos(bamAlign.pos + apath_ref_length(bamAlign.path));
        candidates.push_back(GetSplitSVCandidate(opt,bamRead.target_id(),clipPos,clipPos,isComplex));
    }
}



/// get SV candidates from anomalous read pairs
static
void
getSVCandidatesFromPair(
    const ReadScannerOptions& opt,
    const SVLocusScanner::CachedReadGroupStats& rstats,
    const bam_record& localRead,
    const SimpleAlignment& localAlign,
    const bam_record* remoteReadPtr,
    std::vector<SVCandidate>& candidates)
{
    if (localRead.is_unmapped() || localRead.is_mate_unmapped()) return;

    // update localEvidenceRange:
    const unsigned readSize(apath_read_length(localAlign.path));
    const unsigned localRefLength(apath_ref_length(localAlign.path));

    unsigned thisReadNoninsertSize(0);
    if (localAlign.is_fwd_strand)
    {
        thisReadNoninsertSize=(readSize-apath_read_trail_size(localAlign.path));
    }
    else
    {
        thisReadNoninsertSize=(readSize-apath_read_lead_size(localAlign.path));
    }

    SVCandidate sv;

    SVBreakend& localBreakend(sv.bp1);
    SVBreakend& remoteBreakend(sv.bp2);

    localBreakend.readCount = 1;

    // if remoteRead is not available, estimate mate localRead size to be same as local,
    // and assume no clipping on mate localRead:
    unsigned remoteReadNoninsertSize(readSize);
    unsigned remoteRefLength(readSize);

    if (NULL != remoteReadPtr)
    {
        // if remoteRead is available, we can more accurately determine the size:
        const bam_record& remoteRead(*remoteReadPtr);

        ALIGNPATH::path_t remoteApath;
        bam_cigar_to_apath(remoteRead.raw_cigar(),remoteRead.n_cigar(),remoteApath);

        const unsigned remoteReadSize(apath_read_length(remoteApath));
        remoteRefLength = (apath_ref_length(remoteApath));

        if (remoteRead.is_fwd_strand())
        {
            remoteReadNoninsertSize=(remoteReadSize-apath_read_trail_size(remoteApath));
        }
        else
        {
            remoteReadNoninsertSize=(remoteReadSize-apath_read_lead_size(remoteApath));
        }

        remoteBreakend.readCount = 1;

        localBreakend.pairCount = 1;
        remoteBreakend.pairCount = 1;
    }

    // this is only designed to be valid when reads are on the same chrom with default orientation:
    known_pos_range2 insertRange;

    const pos_t totalNoninsertSize(thisReadNoninsertSize+remoteReadNoninsertSize);
    const pos_t breakendSize(std::max(
                                 static_cast<pos_t>(opt.minPairBreakendSize),
                                 static_cast<pos_t>(rstats.breakendRegion.max-totalNoninsertSize)));

    {
        localBreakend.interval.tid = (localRead.target_id());

        const pos_t startRefPos(localRead.pos()-1);
        const pos_t endRefPos(startRefPos+localRefLength);
        // expected breakpoint range is from the end of the localRead alignment to the (probabilistic) end of the fragment:
        if (localRead.is_fwd_strand())
        {
            localBreakend.state = SVBreakendState::RIGHT_OPEN;
            localBreakend.interval.range.set_begin_pos(endRefPos);
            localBreakend.interval.range.set_end_pos(endRefPos + breakendSize);

            insertRange.set_begin_pos(endRefPos);
        }
        else
        {
            localBreakend.state = SVBreakendState::LEFT_OPEN;
            localBreakend.interval.range.set_end_pos(startRefPos);
            localBreakend.interval.range.set_begin_pos(startRefPos - breakendSize);

            insertRange.set_end_pos(startRefPos);
        }
    }

    // get remote breakend estimate:
    {
        remoteBreakend.interval.tid = (localRead.mate_target_id());

        const pos_t startRefPos(localRead.mate_pos()-1);
        pos_t endRefPos(startRefPos+remoteRefLength);
        if (localRead.is_mate_fwd_strand())
        {
            remoteBreakend.state = SVBreakendState::RIGHT_OPEN;
            remoteBreakend.interval.range.set_begin_pos(endRefPos);
            remoteBreakend.interval.range.set_end_pos(endRefPos + breakendSize);

            insertRange.set_begin_pos(endRefPos);
        }
        else
        {
            remoteBreakend.state = SVBreakendState::LEFT_OPEN;
            remoteBreakend.interval.range.set_end_pos(startRefPos);
            remoteBreakend.interval.range.set_begin_pos(startRefPos - breakendSize);

            insertRange.set_end_pos(startRefPos);
        }
    }

#ifdef DEBUG_SCANNER
    static const std::string logtag("getSVCandidatesFromPair");
    log_os << logtag << " evaluating sv: " << sv << "\n";
#endif


    // check if read pair separation is non-anomalous after accounting for read alignments:
    if (localRead.target_id() == localRead.mate_target_id())
    {
        if (localRead.is_fwd_strand() != localRead.is_mate_fwd_strand())
        {
            // get length of fragment after accounting for any variants described directly in either read alignment:
            const pos_t cigarAdjustedFragmentSize(totalNoninsertSize + (insertRange.end_pos() - insertRange.begin_pos()));

            const bool isLargeFragment(cigarAdjustedFragmentSize > (rstats.properPair.max + opt.minCandidateVariantSize));

            // this is an arbitrary point to start officially tagging 'outties' -- for now  we just want to avoid conventional small fragments from FFPE
            const bool isOuttie(cigarAdjustedFragmentSize < 0);

            if (! (isLargeFragment || isOuttie)) return;
        }
        else
        {
            if (std::abs(localRead.template_size()) <= (rstats.properPair.max + opt.minCandidateVariantSize)) return;
        }
    }

    candidates.push_back(sv);
}



/// get SV candidates from shadow/singleton pairs
/// look for singletons, create candidateSV around conf. interval of shadow position
/// cache singletons? might be needed to remove poor quality shadows.
/// should be able to re-use code, follow soft-clipping example.
#if 0
static
void
getSVCandidatesFromShadow(
    const ReadScannerOptions& opt,
    const SVLocusScanner::CachedReadGroupStats& rstats,
    const bam_record& localRead,
    const SimpleAlignment& localAlign,
    const bam_record* remoteReadPtr,
    std::vector<SVCandidate>& candidates)
{
    if (NULL == remoteReadPtr)
    {
        // we want info on both reads
        return;
    }
    static const bool isComplex(true);
    const bam_record& remoteRead(*remoteReadPtr);
    const SimpleAlignment remoteAlign(remoteRead);

    pos_t singletonGenomePos(0);
    int targetId(0);
    if (localRead.is_mate_unmapped())
    {
        // remote read is shadow candidate
        if (!isGoodShadow(remoteRead,localRead.map_qual(),localRead.qname(),opt.minSingletonMapqGraph))
        {
            return;
        }
        singletonGenomePos = localAlign.pos;
        targetId = remoteRead.target_id();
    }
    else if (localRead.is_unmapped())
    {
        // local is shadow candidate
        if (!isGoodShadow(localRead,remoteRead.map_qual(),remoteRead.qname(),opt.minSingletonMapqGraph))
        {
            return;
        }
        singletonGenomePos = remoteAlign.pos;
        targetId = localRead.target_id();
    }
    else
    {
        // none unmapped, skip this one
        return;
    }
    const pos_t properPairRangeOffset = static_cast<int>(rstats.properPair.min + (rstats.properPair.max-rstats.properPair.min)/2);
    const pos_t shadowGenomePos = singletonGenomePos + properPairRangeOffset;
    candidates.push_back(GetSplitSVCandidate(opt,targetId,shadowGenomePos,shadowGenomePos,isComplex));
}
#endif



/// scan read record (and optionally its mate record) for SV evidence.
//
/// note that estimation is improved by the mate record (because we have the mate cigar string in this case)
///
static
void
getReadBreakendsImpl(
    const ReadScannerOptions& opt,
    const ReadScannerDerivOptions& dopt,
    const SVLocusScanner::CachedReadGroupStats& rstats,
    const bam_record& localRead,
    const bam_record* remoteReadPtr,
    const chromMap_t& chromToIndex,
    std::vector<SVCandidate>& candidates,
    known_pos_range2& localEvidenceRange)
{
    using namespace illumina::common;

    /// TODO: can't handle these yet, but plan to soon:
    //if (localRead.is_mate_unmapped()) return;

    /// get soem basic derived information from the bam_record:
    const ChromAlignment localAlign(localRead);

    // - process any large indels in the localRead:
    getSVCandidatesFromReadIndels(opt, dopt, localAlign, candidates);

#ifdef DEBUG_SCANNER
    static const std::string logtag("getReadBreakendsImpl");
    log_os << logtag << " post-indels candidate_size: " << candidates.size() << "\n";
#endif

    // - process soft-clip in the localRead:
    getSVCandidatesFromReadClip(opt, localRead, localAlign, candidates);

#ifdef DEBUG_SCANNER
    log_os << logtag << " post-clip candidate_size: " << candidates.size() << "\n";
#endif

    // TODO: add semi-aligned read processing
    //
    // CTS: temporarily mark out semi-aligned read input pending review of results with corrected qual offset
    //
    if (false)
    {
        getSVCandidatesFromSemiAligned(opt, localRead, localAlign, candidates);
    }

    /// - process split/SA reads:
    getSACandidatesFromRead(dopt, localRead, localAlign, chromToIndex, candidates);

#ifdef DEBUG_SCANNER
    log_os << logtag << " post-split read candidate_size: " << candidates.size() << "\n";
#endif

    // TODO: process shadow reads
    //getSVCandidatesFromShadow(opt, rstats, localRead, localAlign,remoteReadPtr,candidates);

    // - process anomalous read pairs:
    getSVCandidatesFromPair(opt, rstats, localRead, localAlign, remoteReadPtr, candidates);

#ifdef DEBUG_SCANNER
    log_os << logtag << " post-pair candidate_size: " << candidates.size() << "\n";
#endif

    // update localEvidence range:
    // note this is only used if candidates were added, so there's no harm is setting it every time:
    const unsigned localRefLength(apath_ref_length(localAlign.path));
    const pos_t startRefPos(localRead.pos()-1);
    const pos_t endRefPos(startRefPos+localRefLength);

    localEvidenceRange.set_range(startRefPos,endRefPos);


    /// final chance to QC candidate set:
    ///
    BOOST_FOREACH(const SVCandidate& sv, candidates)
    {
        bool isInvalidTid(false);
        if (sv.bp1.interval.tid < 0)
        {
            isInvalidTid=true;
        }
        else if (sv.bp2.state != SVBreakendState::UNKNOWN)
        {
            if (sv.bp2.interval.tid < 0)
            {
                isInvalidTid=true;
            }
        }

        if (isInvalidTid)
        {
            std::ostringstream oss;
            oss << "SVbreakend has unknown chromosome id in candidate sv.\n"
                << "\tlocal_bam_record: " <<  localRead << "\n"
                << "\tremote_bam record: ";
            if (NULL==remoteReadPtr)
            {
                oss << "NONE";
            }
            else
            {
                oss << (*remoteReadPtr);
            }
            oss << "\n"
                << "\tSVCandidate: " << sv << "\n";
            BOOST_THROW_EXCEPTION(LogicException(oss.str()));

        }
    }
}



/// Create an SVLocus for each potential SV event supported by the BAM record
///
/// the loci count should almost always be one (or, depending on input filtration, zero).
/// multiple suggested loci from one read is more of a theoretical possibility than an
/// expectation.
///
static
void
getSVLociImpl(
    const ReadScannerOptions& opt,
    const ReadScannerDerivOptions& dopt,
    const SVLocusScanner::CachedReadGroupStats& rstats,
    const bam_record& bamRead,
    const chromMap_t& chromToIndex,
    std::vector<SVLocus>& loci)
{
    using namespace illumina::common;

    loci.clear();
    std::vector<SVCandidate> candidates;
    known_pos_range2 localEvidenceRange;

    getReadBreakendsImpl(opt, dopt, rstats, bamRead, NULL, chromToIndex, candidates, localEvidenceRange);

#ifdef DEBUG_SCANNER
    static const std::string logtag("getSVLociImpl");
    log_os << logtag << " candidate_size: " << candidates.size() << "\n";
#endif

    // translate SVCandidate to a simpler form for use
    // in the SV locus graph:
    BOOST_FOREACH(const SVCandidate& cand, candidates)
    {
        const SVBreakend& localBreakend(cand.bp1);
        const SVBreakend& remoteBreakend(cand.bp2);

        const bool isComplex((localBreakend.state == SVBreakendState::COMPLEX) &&
                             (remoteBreakend.state == SVBreakendState::UNKNOWN));

        if ((0==localBreakend.interval.range.size()) ||
            ((! isComplex) && (0==remoteBreakend.interval.range.size())))
        {
            std::ostringstream oss;
            oss << "Unexpected breakend pattern proposed from bam record.\n"
                << "\tlocal_breakend: " << localBreakend << "\n"
                << "\tremote_breakend: " << remoteBreakend << "\n"
                << "\tbam_record: " << bamRead << "\n";
            BOOST_THROW_EXCEPTION(LogicException(oss.str()));
        }

        // determine the evidence weight of this candidate:
        unsigned localEvidenceWeight(0);
        unsigned remoteEvidenceWeight(0);

        if (localBreakend.splitCount != 0)
        {
            localEvidenceWeight = SVObservationWeights::internalReadEvent;
            if (remoteBreakend.splitCount != 0)
            {
                remoteEvidenceWeight = SVObservationWeights::internalReadEvent;
            }
        }
        else if (localBreakend.readCount != 0)
        {
            bool isClose(false);
            if (is_innie_pair(bamRead))
            {
                isClose = (std::abs(bamRead.template_size()) < rstats.minFarFragmentSize);
            }

            unsigned thisWeight(SVObservationWeights::readPair);
            if (isClose) thisWeight = SVObservationWeights::closeReadPair;

            localEvidenceWeight = thisWeight;
            if (remoteBreakend.readCount != 0)
            {
                remoteEvidenceWeight = thisWeight;
            }
        }

        // finally, create the graph locus:
        SVLocus locus;
        // set local breakend estimate:
        const NodeIndexType localBreakendNode(locus.addNode(localBreakend.interval));
        locus.setNodeEvidence(localBreakendNode,localEvidenceRange);

        if (isComplex)
        {
            locus.linkNodes(localBreakendNode,localBreakendNode,localEvidenceWeight);
        }
        else
        {
            // set remote breakend estimate:
            const NodeIndexType remoteBreakendNode(locus.addNode(remoteBreakend.interval));
            locus.linkNodes(localBreakendNode,remoteBreakendNode,localEvidenceWeight,remoteEvidenceWeight);

            locus.mergeSelfOverlap();
        }

        loci.push_back(locus);
    }
}



SVLocusScanner::
SVLocusScanner(
    const ReadScannerOptions& opt,
    const std::string& statsFilename,
    const std::vector<std::string>& alignmentFilename) :
    _opt(opt),
    _dopt(opt)
{
    using namespace illumina::common;

    // pull in insert stats:
    _rss.load(statsFilename.c_str());

    // cache the insert stats we'll be looking up most often:
    BOOST_FOREACH(const std::string& alignFile, alignmentFilename)
    {
        const boost::optional<unsigned> index(_rss.getGroupIndex(alignFile));
        if (! index)
        {
            std::ostringstream oss;
            oss << "Can't find alignment file in pre-computed fragment stats.\n"
                << "\tstatsFile: " << statsFilename << "\n"
                << "\talignFile: " << alignFile << "\n";
            BOOST_THROW_EXCEPTION(LogicException(oss.str()));
        }
        const ReadGroupStats rgs(_rss.getStats(*index));

        _stats.resize(_stats.size()+1);
        CachedReadGroupStats& stat(_stats.back());
        {
            Range& breakend(stat.breakendRegion);
            breakend.min=rgs.fragStats.quantile(_opt.breakendEdgeTrimProb);
            breakend.max=rgs.fragStats.quantile((1-_opt.breakendEdgeTrimProb));

            if (breakend.min<0.) breakend.min = 0;
            assert(breakend.max>0.);
        }
        {
            Range& ppair(stat.properPair);
            ppair.min=rgs.fragStats.quantile(_opt.properPairTrimProb);
            ppair.max=rgs.fragStats.quantile((1-_opt.properPairTrimProb));

            if (ppair.min<0.) ppair.min = 0;

            assert(ppair.max>0.);
        }
        {
            Range& evidence(stat.evidencePair);
            evidence.min=rgs.fragStats.quantile(_opt.evidenceTrimProb);
            evidence.max=rgs.fragStats.quantile((1-_opt.evidenceTrimProb));

            if (evidence.min<0.) evidence.min = 0;

            assert(evidence.max>0.);
        }

        stat.minFarFragmentSize = static_cast<int>(stat.properPair.max*SVObservationWeights::closePairFactor);
    }
}



bool
SVLocusScanner::
isReadFiltered(const bam_record& bamRead) const
{
    if (bamRead.is_filter()) return true;
    if (bamRead.is_dup()) return true;
    if (bamRead.is_secondary()) return true;
    if (bamRead.map_qual() < _opt.minMapq) return true;
    return false;
}



bool
SVLocusScanner::
isProperPair(
    const bam_record& bamRead,
    const unsigned defaultReadGroupIndex) const
{
    if (! is_innie_pair(bamRead)) return false;

    const Range& ppr(_stats[defaultReadGroupIndex].properPair);
    const int32_t fragmentSize(std::abs(bamRead.template_size()));

    /// we're seeing way to much large fragment garbage in cancers to use the normal proper pair criteria, push the max fragment size out a bit for now:
    static const float maxAnomFactor(1.5);
    if ((fragmentSize > static_cast<int32_t>(maxAnomFactor*ppr.max)) || (fragmentSize < ppr.min)) return false;

    return true;
}



bool
SVLocusScanner::
isLargeFragment(
    const bam_record& bamRead,
    const unsigned defaultReadGroupIndex) const
{
    if (bamRead.target_id() != bamRead.mate_target_id()) return true;

    const Range& ppr(_stats[defaultReadGroupIndex].properPair);
    const int32_t fragmentSize(std::abs(bamRead.template_size()));
    return (fragmentSize > ppr.max);
}


bool
SVLocusScanner::
isLocalAssemblyEvidence(
    const bam_record& bamRead) const
{
    using namespace ALIGNPATH;

    {
        // TODO: (1) double check semi-aligned thresholds with fixed qual offsets
        // TODO: For the semi-aligned test to be effective, we need to convert the cigar string to contain match/mis-match infomration first
        //  if (isSemiAligned(bamRead,_opt.minSemiAlignedScoreGraph)) return true;
    }

    const SimpleAlignment bamAlign(bamRead);
    //
    // large indel already in cigar string
    //
    BOOST_FOREACH(const path_segment& ps, bamAlign.path)
    {
        if (ps.type == INSERT || ps.type == DELETE)
        {
            if (ps.length>=_opt.minCandidateVariantSize) return true;
        }
    }

    //
    // soft-clipping:
    //
    {

        unsigned leadingClipLen(0), trailingClipLen(0);
        getSVBreakendCandidateClip(bamRead, bamAlign.path, leadingClipLen, trailingClipLen);

        if ((leadingClipLen >= _opt.minSoftClipLen) || (trailingClipLen >= _opt.minSoftClipLen)) return true;
    }

    return false;
}




void
SVLocusScanner::
getSVLoci(
    const bam_record& bamRead,
    const unsigned defaultReadGroupIndex,
    const std::map<std::string, int32_t>& chromToIndex,
    std::vector<SVLocus>& loci) const
{
    loci.clear();

    const CachedReadGroupStats& rstats(_stats[defaultReadGroupIndex]);
    getSVLociImpl(_opt, _dopt, rstats, bamRead, chromToIndex, loci);
}



void
SVLocusScanner::
getBreakendPair(
    const bam_record& localRead,
    const bam_record* remoteReadPtr,
    const unsigned defaultReadGroupIndex,
    const  std::map<std::string, int32_t>& chromToIndex,
    std::vector<SVCandidate>& candidates) const
{
    const CachedReadGroupStats& rstats(_stats[defaultReadGroupIndex]);

    // throw evidence range away in this case
    known_pos_range2 evidenceRange;
    getReadBreakendsImpl(_opt, _dopt, rstats, localRead, remoteReadPtr, chromToIndex, candidates, evidenceRange);
}
