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
/// \author Chris Saunders
/// \author Ole Schulz-Trieglaff
/// \author Bret Barnes
///

#include "blt_util/align_path_bam_util.hh"
#include "blt_util/align_path_util.hh"
#include "blt_util/bam_record_util.hh"

#include "blt_util/parse_util.hh"
#include "blt_util/string_util.hh"
#include "common/Exceptions.hh"
#include "manta/SVCandidateUtil.hh"
#include "manta/SVLocusScanner.hh"
#include "manta/SVLocusScannerSemiAligned.hh"

#include "boost/foreach.hpp"


// #define DEBUG_SCANNER

//#define DEBUG_IS_SHADOW

#ifdef DEBUG_SCANNER
#include "blt_util/log.hh"
#include <iostream>
#endif


/// used for classifying fragments based on size so that they can be treated differently
///
namespace FragmentSizeType
{
static const float closePairFactor(4); ///< fragments within this factor of the minimum size cutoff are treated as 'close' pairs and receive a modified evidence count
static const float veryClosePairFactor(1.5); ///< fragments within this factor of the minimum size cutoff are treated as 'reallyClose' pairs and receive a modified evidence count
static const float maxNormalFactor(1.5);

static const float minLargeEventRegionFactor(10);
static const float maxLargeEventRegionFactor(20);


static
index_t
classifySize(
    const SVLocusScanner::CachedReadGroupStats& rgStats,
    const int fragmentSize)
{
    if (fragmentSize < rgStats.properPair.min) return COMPRESSED;
    if (fragmentSize > rgStats.properPair.max)
    {
        if (fragmentSize < rgStats.minDistantFragmentSize) return CLOSE;
        return DISTANT;
    }
    return NORMAL;
}

static
bool
isLarge(const index_t i)
{
    switch (i)
    {
    case NORMAL:
    case COMPRESSED:
        return false;
    default:
        return true;
    }
}
}



/// temporary object used to track the insertion
/// of SVCandidates into the pool:
///
struct TrackedCandidates
{
    TrackedCandidates(
        std::vector<SVObservation>& candidates,
        TruthTracker& truthTracker) :
        data(candidates),
        _truthTracker(truthTracker)
    {}

    void
    push_back(const SVObservation& obs)
    {
        data.push_back(obs);
        _truthTracker.addObservation(data.back());
    }

    unsigned
    size() const
    {
        return data.size();
    }

    std::vector<SVObservation>& data;
private:
    TruthTracker& _truthTracker;
};




static
SVObservation
GetSplitSVCandidate(
    const ReadScannerDerivOptions& dopt,
    const int32_t alignTid,
    const pos_t leftPos,
    const pos_t rightPos,
    const SVEvidenceType::index_t& svSource,
    const bool isComplex = false)
{
    SVObservation sv;
    SVBreakend& localBreakend(sv.bp1);
    SVBreakend& remoteBreakend(sv.bp2);

    localBreakend.interval.tid = alignTid;
    remoteBreakend.interval.tid = alignTid;

    localBreakend.lowresEvidence.add(svSource);
    sv.evtype = svSource;

    if (! isComplex)
    {
        remoteBreakend.lowresEvidence.add(svSource);
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
    const SimpleAlignment& align,
    SVBreakend& breakend)
{
    // Need to use the match descriptors to determine if
    // we are upstream clipped or downstream clipped.
    // Below is the logic to convert these  to breakend candidates:
    //
    // Upstream => LEFT_OPEN
    // DownStream => RIGHT_OPEN
    //

    const bool isUpstream(isCigarUpstream(align.path));

    if (isUpstream)
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



/// get SV candidates from SA-tag split-read alignment
static
SVObservation
GetSplitSACandidate(
    const ReadScannerDerivOptions& dopt,
    const bam_record& localRead,
    const SimpleAlignment& localAlign,
    const SimpleAlignment& remoteAlign)
{
    using namespace SVEvidenceType;
    static const index_t svSource(SPLIT_ALIGN);

    SVObservation sv;

    sv.bp1.lowresEvidence.add(svSource);
    sv.bp2.lowresEvidence.add(svSource);
    sv.evtype = svSource;

    updateSABreakend(dopt, localAlign, sv.bp1);
    updateSABreakend(dopt, remoteAlign, sv.bp2);

    bool remoteIsSecond = !ALIGNPATH::is_clipped_front(localAlign.path); // todo Need smarter check which partial alignment is first
    bool isReadFw = (localRead.is_first() == localRead.is_fwd_strand());
    if (isReadFw == remoteIsSecond)
    {
        sv.fwReads += 1;
    }
    else
    {
        sv.rvReads += 1;
    }

    return sv;
}



typedef std::map<std::string, int32_t> chromMap_t;



static
void
getSACandidatesFromRead(
    const ReadScannerOptions& opt,
    const ReadScannerDerivOptions& dopt,
    const bam_record& localRead,
    const SimpleAlignment& localAlign,
    const chromMap_t& chromToIndex,
    TrackedCandidates& candidates)
{
    using namespace ALIGNPATH;

    std::vector<std::string> saVec;
    {
        static const char satag[] = {'S','A'};
        const char* saStr(localRead.get_string_tag(satag));
        if (NULL == saStr) return;

        split_string(saStr, ';', saVec);
        if ( (! saVec.empty()) && saVec.back().empty())
        {
            saVec.pop_back();
        }
    }

    // Only handle a single split alignment right now.
    // In the future we could sort the SA tags by order on the template, possibly
    // also removing segments that map to two different areas,

    if (saVec.size() > 1) return;

    SimpleAlignment remoteAlign;

    BOOST_FOREACH(const std::string& sa, saVec)
    {
#ifdef DEBUG_SCANNER
        log_os << "SA STRING: " << sa << "\n";
#endif
        std::vector<std::string> saDat;
        split_string(sa, ',', saDat);

        assert((saDat.size() == 6) && "Unexpected number of SA tag values");

        /// filter split reads with low MAPQ:
        const unsigned saMapq(illumina::blt_util::parse_unsigned_str(saDat[4]));
        if (saMapq < opt.minMapq) continue;

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

        candidates.push_back(GetSplitSACandidate(dopt, localRead, localAlign, remoteAlign));
    }
}



static
void
getSVCandidatesFromReadIndels(
    const ReadScannerOptions& opt,
    const ReadScannerDerivOptions& dopt,
    const SimpleAlignment& align,
    TrackedCandidates& candidates)
{
    using namespace SVEvidenceType;
    static const index_t svSource(CIGAR);

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
                    candidates.push_back(GetSplitSVCandidate(dopt, align.tid, refHeadPos, refHeadPos, svSource, isComplex));
                }
            }
        }
        else if (isSwapStart)
        {
            const swap_info sinfo(align.path,pathIndex);
            if ((sinfo.delete_length >= opt.minCandidateVariantSize) || (sinfo.insert_length >= opt.minCandidateVariantSize))
            {
                candidates.push_back(GetSplitSVCandidate(dopt, align.tid, refHeadPos, refHeadPos+sinfo.delete_length, svSource));
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
                    candidates.push_back(GetSplitSVCandidate(dopt, align.tid, refHeadPos, refHeadPos+ps.length, svSource));
                }
            }
            else if (ps.type == INSERT)
            {
                if (ps.length >= opt.minCandidateVariantSize)
                {
                    candidates.push_back(GetSplitSVCandidate(dopt, align.tid, refHeadPos, refHeadPos, svSource));
                }
            }
        }

        for (unsigned i(0); i<nPathSegments; ++i)
        {
            increment_path(align.path, pathIndex, readOffset, refHeadPos);
        }
    }
}



static
void
getSVCandidatesFromSemiAligned(
    const ReadScannerOptions& opt,
    const bam_record& bamRead,
    const SimpleAlignment& bamAlign,
    const reference_contig_segment& refSeq,
    TrackedCandidates& candidates)
{
    unsigned leadingMismatchLen(0);
    unsigned trailingMismatchLen(0);
    pos_t leadingRefPos(0), trailingRefPos(0);
    getSVBreakendCandidateSemiAligned(bamRead, bamAlign, refSeq,
                                      leadingMismatchLen, leadingRefPos,
                                      trailingMismatchLen, trailingRefPos);

    if ((leadingMismatchLen + trailingMismatchLen) >= bamRead.read_size()) return;

    using namespace SVEvidenceType;
    static const index_t svSource(SEMIALIGN);

    // semi-aligned reads don't define a full hypothesis, so they're always evidence for a 'complex' ie. undefined, event
    // in a fashion analogous to clipped reads
    static const bool isComplex(true);

    if (leadingMismatchLen >= opt.minSemiAlignedMismatchLen)
    {
        const pos_t pos(leadingRefPos);
        candidates.push_back(GetSplitSVCandidate(opt,bamRead.target_id(),pos,pos,svSource,isComplex));
    }

    if (trailingMismatchLen >= opt.minSemiAlignedMismatchLen)
    {
        const pos_t pos(trailingRefPos);
        candidates.push_back(GetSplitSVCandidate(opt,bamRead.target_id(),pos,pos,svSource,isComplex));
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
    TrackedCandidates& candidates)
{
    using namespace SVEvidenceType;
    static const index_t svLocalPair(LOCAL_PAIR);
    static const index_t svPair(PAIR);

    if (! localRead.is_paired()) return;

    if (localRead.is_unmapped() || localRead.is_mate_unmapped()) return;

    // special case typically used for RNA-Seq analysis:
    if (opt.isIgnoreAnomProperPair && localRead.is_proper_pair()) return;

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

    SVObservation sv;
    sv.fwReads = 1; // bp1 is read1, bp2 is read2
    SVBreakend& localBreakend(sv.bp1);
    SVBreakend& remoteBreakend(sv.bp2);

    localBreakend.lowresEvidence.add(svLocalPair);
    sv.evtype = svLocalPair;

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

        remoteBreakend.lowresEvidence.add(svLocalPair);

        localBreakend.lowresEvidence.add(svPair);
        remoteBreakend.lowresEvidence.add(svPair);
        sv.evtype = svPair;
    }

    const pos_t localStartRefPos(localRead.pos()-1);
    const pos_t localEndRefPos(localStartRefPos+localRefLength);
    const pos_t remoteStartRefPos(localRead.mate_pos()-1);
    const pos_t remoteEndRefPos(remoteStartRefPos+remoteRefLength);

    // this is only designed to be valid when reads are on the same chrom with default relative orientation:
    known_pos_range2 insertRange;
    if (localRead.is_fwd_strand())
    {
        insertRange.set_range(localEndRefPos,remoteStartRefPos);
    }
    else
    {
        insertRange.set_range(remoteEndRefPos,localStartRefPos);
    }

    // total the reference span or read1 and read2 (ie. the elements of the
    // fragment that are not part of the insert between the reads)
    //
    const pos_t totalNoninsertSize(thisReadNoninsertSize+remoteReadNoninsertSize);

    // different breakend sizes are used for long-range pairings vs short-ish range deletions,
    // because of occasional long-fragment noise. This ramps from 0 to 1 as we go from short to
    // long deletions sizes:
    double largeEventRegionScale(1.0);

    // check if fragment size is still anomalous after accounting for read alignment patterns:
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

            if (! isOuttie)
            {
                largeEventRegionScale = rstats.largeEventRegionScaler.getScale(cigarAdjustedFragmentSize);
            }
        }
    }


    // set state and interval for each breakend:
    {
        const double breakendRegionMax(
            (largeEventRegionScale*rstats.largeScaleEventBreakendRegion.max) +
            ((1.-largeEventRegionScale)*rstats.breakendRegion.max));

        const pos_t breakendSize(std::max(
                                     static_cast<pos_t>(opt.minPairBreakendSize),
                                     static_cast<pos_t>(breakendRegionMax-totalNoninsertSize)));

        localBreakend.interval.tid = (localRead.target_id());
        // expected breakpoint range is from the end of the localRead alignment to the (probabilistic) end of the fragment:
        if (localRead.is_fwd_strand())
        {
            localBreakend.state = SVBreakendState::RIGHT_OPEN;
            localBreakend.interval.range.set_begin_pos(localEndRefPos);
            localBreakend.interval.range.set_end_pos(localEndRefPos + breakendSize);
        }
        else
        {
            localBreakend.state = SVBreakendState::LEFT_OPEN;
            localBreakend.interval.range.set_end_pos(localStartRefPos);
            localBreakend.interval.range.set_begin_pos(localStartRefPos - breakendSize);
        }

        remoteBreakend.interval.tid = (localRead.mate_target_id());

        if (localRead.is_mate_fwd_strand())
        {
            remoteBreakend.state = SVBreakendState::RIGHT_OPEN;
            remoteBreakend.interval.range.set_begin_pos(remoteEndRefPos);
            remoteBreakend.interval.range.set_end_pos(remoteEndRefPos + breakendSize);
        }
        else
        {
            remoteBreakend.state = SVBreakendState::LEFT_OPEN;
            remoteBreakend.interval.range.set_end_pos(remoteStartRefPos);
            remoteBreakend.interval.range.set_begin_pos(remoteStartRefPos - breakendSize);
        }
    }

#ifdef DEBUG_SCANNER
    static const std::string logtag("getSVCandidatesFromPair");
    log_os << logtag << " evaluating pair sv for inclusion: " << sv << "\n";
#endif

    candidates.push_back(sv);
}


#if 0
/// get SV candidates from shadow/singleton pairs
/// look for singletons, create candidateSV around conf. interval of shadow position
/// cache singletons? might be needed to remove poor quality shadows.
/// should be able to re-use code, follow soft-clipping example.
static
void
getSVCandidatesFromShadow(
    const ReadScannerOptions& opt,
    const SVLocusScanner::CachedReadGroupStats& rstats,
    const bam_record& localRead,
    const SimpleAlignment& localAlign,
    const bam_record* remoteReadPtr,
    TrackedCandidates& candidates)
{
    using namespace SVEvidenceType;
    static const index_t svSource(SHADOW);

    static const bool isComplex(true);
    pos_t singletonGenomePos(0);
    int targetId(0);
    if (NULL == remoteReadPtr)
    {
        if (!localRead.is_unmapped()) return;
        // need to take care of this case
        // need to rely on cached mapq and qname
        return;
        if (!isGoodShadow(localRead,lastMapq,lastQname,opt.minSingletonMapqGraph))
        {
            return;
        }
        singletonGenomePos = localAlign.pos;
        targetId           = localRead.target_id();
    }
    else
    {
        // have both reads, straightforward from here
        const bam_record& remoteRead(*remoteReadPtr);
        const SimpleAlignment remoteAlign(remoteRead);

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
    }
    const pos_t properPairRangeOffset = static_cast<pos_t>(rstats.properPair.min + (rstats.properPair.max-rstats.properPair.min)/2);
    const pos_t shadowGenomePos = singletonGenomePos + properPairRangeOffset;
    candidates.push_back(GetSplitSVCandidate(opt,targetId,shadowGenomePos,shadowGenomePos, svSource, isComplex));
}
#endif



static
void
getSingleReadSVCandidates(
    const ReadScannerOptions& opt,
    const ReadScannerDerivOptions& dopt,
    const bam_record& localRead,
    const SimpleAlignment& localAlign,
    const chromMap_t& chromToIndex,
    const reference_contig_segment& refSeq,
    TrackedCandidates& candidates)
{
    using namespace illumina::common;

    /// TODO: can't handle these yet, but plan to soon:
    //if (localRead.is_mate_unmapped()) return;

    // - process any large indels in the localRead:
    getSVCandidatesFromReadIndels(opt, dopt, localAlign, candidates);
#ifdef DEBUG_SCANNER
    static const std::string logtag("getSingleReadSVCandidates");
    log_os << logtag << " post-indels candidate_size: " << candidates.size() << "\n";
#endif

    // this detects semi-aligned AND soft-clip now:
    getSVCandidatesFromSemiAligned(opt, localRead, localAlign, refSeq,
                                   candidates);
#ifdef DEBUG_SCANNER
    log_os << logtag << " post-semialigned candidate_size: " << candidates.size() << "\n";
#endif

    /// - process split/SA reads:
    getSACandidatesFromRead(opt, dopt, localRead, localAlign, chromToIndex,
                            candidates);
#ifdef DEBUG_SCANNER
    log_os << logtag << " post-split read candidate_size: " << candidates.size() << "\n";
#endif
}



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
    const bam_header_info& bamHeader,
    const reference_contig_segment& localRefSeq,
    const reference_contig_segment* remoteRefSeqPtr,
    TrackedCandidates& candidates,
    known_pos_range2& localEvidenceRange)
{
    using namespace illumina::common;

    const chromMap_t& chromToIndex(bamHeader.chrom_to_index);

    candidates.data.clear();

    /// TODO: can't handle these yet, but plan to soon:
    //if (localRead.is_mate_unmapped()) return;

    /// get some basic derived information from the bam_record:
    const SimpleAlignment localAlign(localRead);

    getSingleReadSVCandidates(opt, dopt, localRead, localAlign, chromToIndex,
                              localRefSeq, candidates);

    if (NULL != remoteReadPtr)
    {
        // run the same check on the read's mate if we have access to it
        assert(NULL != remoteRefSeqPtr);
        const bam_record& remoteRead(*remoteReadPtr);
        const SimpleAlignment remoteAlign(remoteRead);

        getSingleReadSVCandidates(opt, dopt, remoteRead, remoteAlign,
                                  chromToIndex, (*remoteRefSeqPtr),
                                  candidates);
    }

    // process shadows:
    //getSVCandidatesFromShadow(opt, rstats, localRead, localAlign,remoteReadPtr,candidates);

    // - process anomalous read pairs:
    getSVCandidatesFromPair(opt, rstats, localRead, localAlign, remoteReadPtr,
                            candidates);

#ifdef DEBUG_SCANNER
    static const std::string logtag("getReadBreakendsImpl: ");
    log_os << logtag << "post-pair candidate_size: " << candidates.size() << "\n";
#endif

    // update localEvidence range:
    // note this is only used if candidates were added, so there's no harm in setting it every time:
    const unsigned localRefLength(apath_ref_length(localAlign.path));
    const pos_t startRefPos(localRead.pos()-1);
    const pos_t endRefPos(startRefPos+localRefLength);

    localEvidenceRange.set_range(startRefPos,endRefPos);

    const int maxTid(chromToIndex.size());

    /// final chance to QC candidate set:
    ///
    BOOST_FOREACH(const SVCandidate& sv, candidates.data)
    {
        bool isInvalidTid(false);
        if ((sv.bp1.interval.tid < 0) || (sv.bp1.interval.tid >= maxTid))
        {
            isInvalidTid=true;
        }
        else if (sv.bp2.state != SVBreakendState::UNKNOWN)
        {
            if ((sv.bp2.interval.tid < 0) || (sv.bp2.interval.tid >= maxTid))
            {
                isInvalidTid=true;
            }
        }

        bool isInvalidPos(false);
        if (! isInvalidTid)
        {
            // note in the 'off-chromosome edge' test below we check for cases which are obviously way off
            // the edge, but allow for a bit of over-edge mistakes to occur for the circular chromosomes
            //
            static const int offEdgePad(500);
            const pos_t tid1Length(bamHeader.chrom_data[sv.bp1.interval.tid].length);
            if ((sv.bp1.interval.range.end_pos() <= -offEdgePad) || (sv.bp1.interval.range.begin_pos() >= (tid1Length+offEdgePad)))
            {
                isInvalidPos=true;
            }
            else if (sv.bp2.state != SVBreakendState::UNKNOWN)
            {
                const pos_t tid2Length(bamHeader.chrom_data[sv.bp2.interval.tid].length);
                if ((sv.bp2.interval.range.end_pos() <= -offEdgePad) || (sv.bp2.interval.range.begin_pos() >= (tid2Length+offEdgePad)))
                {
                    isInvalidPos=true;
                }
            }
        }

        if (isInvalidTid || isInvalidPos)
        {
            std::ostringstream oss;
            if (isInvalidTid)
            {
                oss << "SVbreakend has unknown or invalid chromosome id in candidate sv.\n";
            }
            else
            {
                oss << "Cannot interpret BAM record: candidate SV breakend from BAM record is off chromosome edge.\n";
            }

            oss << "\tlocal_bam_record: " <<  localRead << "\n"
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
    const bam_header_info& bamHeader,
    const reference_contig_segment& refSeq,
    std::vector<SVLocus>& loci,
    TruthTracker& truthTracker)
{
    using namespace illumina::common;

    loci.clear();
    std::vector<SVObservation> candidates;
    TrackedCandidates trackedCandidates(candidates,truthTracker);
    known_pos_range2 localEvidenceRange;

    getReadBreakendsImpl(opt, dopt, rstats, bamRead, NULL, bamHeader,
                         refSeq, NULL, trackedCandidates, localEvidenceRange);

#ifdef DEBUG_SCANNER
    static const std::string logtag("getSVLociImpl");
    log_os << logtag << " candidate_size: " << candidates.size() << "\n";
#endif

    // translate SVCandidate to a simpler form for use
    // in the SV locus graph:
    BOOST_FOREACH(const SVCandidate& cand, candidates)
    {
        const bool isCandComplex(isComplexSV(cand));

        const SVBreakend& localBreakend(cand.bp1);
        const SVBreakend& remoteBreakend(cand.bp2);

        if ((0==localBreakend.interval.range.size()) ||
            ((! isCandComplex) && (0==remoteBreakend.interval.range.size())))
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

        if (localBreakend.getAnyNonPairCount() != 0)
        {
            localEvidenceWeight = SVObservationWeights::internalReadEvent;
            if (remoteBreakend.getAnyNonPairCount() != 0)
            {
                remoteEvidenceWeight = SVObservationWeights::internalReadEvent;
            }
        }
        else if (localBreakend.getLocalPairCount() != 0)
        {
            bool isClose(false);
            if (is_innie_pair(bamRead))
            {
                isClose = (std::abs(bamRead.template_size()) < rstats.minDistantFragmentSize);
            }

            unsigned thisWeight(SVObservationWeights::readPair);
            if (isClose) thisWeight = SVObservationWeights::closeReadPair;

            localEvidenceWeight = thisWeight;
            if (remoteBreakend.getLocalPairCount() != 0)
            {
                remoteEvidenceWeight = thisWeight;
            }
        }

        // finally, create the graph locus:
        SVLocus locus;
        // set local breakend estimate:
        const NodeIndexType localBreakendNode(locus.addNode(localBreakend.interval));
        locus.setNodeEvidence(localBreakendNode,localEvidenceRange);

        if (isCandComplex)
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



/// compute one of the scanner's fragment ranges:
static
void
setRGRange(
    const SizeDistribution& fragStats,
    const float qprob,
    SVLocusScanner::Range& range)
{
    range.min=fragStats.quantile(qprob);
    range.max=fragStats.quantile((1-qprob));
    if (range.min<0.) range.min = 0;
    assert(range.max>0.);
}



SVLocusScanner::
SVLocusScanner(
    const ReadScannerOptions& opt,
    const std::string& statsFilename,
    const std::vector<std::string>& /*alignmentFilename*/) :
    _opt(opt),
    _dopt(opt)
{
    using namespace illumina::common;

    // pull in insert stats:
    _rss.load(statsFilename.c_str());

    // precompute frequently used insert stats for each rg:
    const unsigned rgCount(_rss.size());
    for (unsigned rgIndex(0); rgIndex<rgCount; rgIndex++)
    {
        /// TODO: add check that the filenames in the stats file are a complete match to alignmentFilename

        const SizeDistribution& rgDistro(getFragSizeDistro(rgIndex));

        _stats.resize(_stats.size()+1);
        CachedReadGroupStats& rgStats(_stats.back());
        setRGRange(rgDistro, _opt.breakendEdgeTrimProb, rgStats.breakendRegion);
        setRGRange(rgDistro, _opt.largeScaleEventBreakendEdgeTrimProb, rgStats.largeScaleEventBreakendRegion);
        setRGRange(rgDistro, _opt.properPairTrimProb, rgStats.properPair);
        setRGRange(rgDistro, _opt.evidenceTrimProb, rgStats.evidencePair);

        rgStats.shadowSearchRange = rgDistro.quantile(1-(_opt.shadowSearchRangeProb))*_opt.shadowSearchRangeFactor;

        assert(rgStats.shadowSearchRange > 0);

        rgStats.minVeryCloseFragmentSize = static_cast<int>(rgStats.properPair.max*FragmentSizeType::maxNormalFactor);
        rgStats.minCloseFragmentSize = static_cast<int>(rgStats.properPair.max*FragmentSizeType::veryClosePairFactor);
        rgStats.minDistantFragmentSize = static_cast<int>(rgStats.properPair.max*FragmentSizeType::closePairFactor);

        assert(rgStats.minDistantFragmentSize > rgStats.properPair.max);

        //rgStats.veryCloseEventScaler.init(rgStats.minVeryCloseFragmentSize, rgStats.minCloseFragmentSize);

        const int largeEventRegionMin(rgStats.properPair.max*FragmentSizeType::minLargeEventRegionFactor);
        const int largeEventRegionMax(rgStats.properPair.max*FragmentSizeType::maxLargeEventRegionFactor);

        rgStats.largeEventRegionScaler.init(largeEventRegionMin, largeEventRegionMax);
    }
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



FragmentSizeType::index_t
SVLocusScanner::
getFragmentSizeType(
    const bam_record& bamRead,
    const unsigned defaultReadGroupIndex) const
{
    using namespace FragmentSizeType;
    if (bamRead.target_id() != bamRead.mate_target_id()) return DISTANT;
    const int32_t fragmentSize(std::abs(bamRead.template_size()));
    return classifySize(_stats[defaultReadGroupIndex], fragmentSize);
}


bool
SVLocusScanner::
isLargeFragment(
    const bam_record& bamRead,
    const unsigned defaultReadGroupIndex) const
{
    if (! bamRead.is_paired()) return false;
    return FragmentSizeType::isLarge(getFragmentSizeType(bamRead,defaultReadGroupIndex));
}



bool
SVLocusScanner::
isNonCompressedAnomalous(
    const bam_record& bamRead,
    const unsigned defaultReadGroupIndex) const
{
    const bool isAnomalous(! isProperPair(bamRead,defaultReadGroupIndex));
    const bool isInnie(is_innie_pair(bamRead));
    const bool isLarge(isLargeFragment(bamRead,defaultReadGroupIndex));

    // exclude innie read pairs which are anomalously short:
    return (isAnomalous && ((! isInnie) || isLarge));
}



bool
SVLocusScanner::
isLocalAssemblyEvidence(
    const bam_record& bamRead,
    const reference_contig_segment& refSeq) const
{
    using namespace ALIGNPATH;

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
    // semi-aligned AND soft-clipped read ends:
    //
    {
        unsigned leadingMismatchLen(0), trailingMismatchLen(0);
        getSVBreakendCandidateSemiAligned(bamRead, bamAlign, refSeq, leadingMismatchLen, trailingMismatchLen);
        if ((leadingMismatchLen >= _opt.minSemiAlignedMismatchLen) || (trailingMismatchLen >= _opt.minSemiAlignedMismatchLen))
        {
            return true;
        }
    }

    return false;
}




void
SVLocusScanner::
getSVLoci(
    const bam_record& bamRead,
    const unsigned defaultReadGroupIndex,
    const bam_header_info& bamHeader,
    const reference_contig_segment& refSeq,
    std::vector<SVLocus>& loci,
    TruthTracker& truthTracker) const
{
    loci.clear();

    const CachedReadGroupStats& rstats(_stats[defaultReadGroupIndex]);
    getSVLociImpl(_opt, _dopt, rstats, bamRead, bamHeader, refSeq, loci,
                  truthTracker);
    //lastQname = bamRead.qname();
    //lastMapq  = bamRead.map_qual();
}



void
SVLocusScanner::
getBreakendPair(
    const bam_record& localRead,
    const bam_record* remoteReadPtr,
    const unsigned defaultReadGroupIndex,
    const bam_header_info& bamHeader,
    const reference_contig_segment& localRefSeq,
    const reference_contig_segment* remoteRefSeqPtr,
    std::vector<SVObservation>& candidates,
    TruthTracker& truthTracker) const
{
    const CachedReadGroupStats& rstats(_stats[defaultReadGroupIndex]);

    // throw evidence range away in this case
    TrackedCandidates trackedCandidates(candidates, truthTracker);
    known_pos_range2 evidenceRange;
    getReadBreakendsImpl(_opt, _dopt, rstats, localRead, remoteReadPtr,
                         bamHeader, localRefSeq, remoteRefSeqPtr,
                         trackedCandidates, evidenceRange);
}
