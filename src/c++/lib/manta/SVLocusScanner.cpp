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
// <https://github.com/downloads/sequencing/licenses/>.
//

///
/// \author Chris Saunders
///

#include "manta/SVLocusScanner.hh"
#include "blt_util/align_path_bam_util.hh"
#include "blt_util/align_path_util.hh"
#include "blt_util/log.hh"
#include "common/Exceptions.hh"

#include "boost/foreach.hpp"



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



static
SVCandidate
GetSplitSVCandidate(
    const ReadScannerOptions& opt,
    const int32_t alignTid,
    const pos_t leftPos,
    const pos_t rightPos)
{
    SVCandidate sv;
    SVBreakend& localBreakend(sv.bp1);
    SVBreakend& remoteBreakend(sv.bp2);

    localBreakend.splitCount++;
    remoteBreakend.splitCount++;

    localBreakend.state = SVBreakendState::RIGHT_OPEN;
    remoteBreakend.state = SVBreakendState::LEFT_OPEN;

    localBreakend.interval.tid = alignTid;
    remoteBreakend.interval.tid = alignTid;

    const pos_t beforeBreakend(opt.minPairBreakendSize/2);
    const pos_t afterBreakend(opt.minPairBreakendSize-beforeBreakend);

    localBreakend.interval.range.set_begin_pos(std::max(0,leftPos-beforeBreakend));
    localBreakend.interval.range.set_end_pos(leftPos+afterBreakend);

    remoteBreakend.interval.range.set_begin_pos(std::max(0,rightPos-beforeBreakend));
    remoteBreakend.interval.range.set_end_pos(rightPos+afterBreakend);

    return sv;
}



static
void
getSVCandidatesFromRead(
    const ReadScannerOptions& opt,
    const SimpleAlignment& align,
    const int32_t alignTid,
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

        assert(ps.type != SKIP);
        assert(! (isEdgeSegment && isSwapStart));

        unsigned nPathSegments(1); // number of path segments consumed
        if (isEdgeSegment)
        {
            // edge inserts are allowed for intron adjacent and grouper reads, edge deletions for intron adjacent only

            if (ps.type == INSERT)
            {
                // ignore for now...
            }
            else if (ps.type == SOFT_CLIP)
            {
                // ignore for now...
            }
        }
        else if (isSwapStart)
        {
            const swap_info sinfo(align.path,pathIndex);
            if (sinfo.delete_length >= opt.minCandidateIndelSize)
            {
                candidates.push_back(GetSplitSVCandidate(opt,alignTid,refHeadPos,refHeadPos+sinfo.delete_length));
            }

            nPathSegments = sinfo.n_seg;
        }
        else if (is_segment_type_indel(align.path[pathIndex].type))
        {
            // regular indel:

            if (ps.type == DELETE)
            {
                if (align.path[pathIndex].length >= opt.minCandidateIndelSize)
                {
                    candidates.push_back(GetSplitSVCandidate(opt,alignTid,refHeadPos,refHeadPos+align.path[pathIndex].length));
                }
            }

            // ignore other indel types for now...
        }

        for (unsigned i(0); i<nPathSegments; ++i)
        {
            increment_path(align.path,pathIndex,readOffset,refHeadPos);
        }
    }
}



static
void
getReadBreakendsImpl(
    const ReadScannerOptions& opt,
    const SVLocusScanner::CachedReadGroupStats& rstats,
    const bam_record& localRead,
    const bam_record* remoteReadPtr,
    std::vector<SVCandidate>& candidates,
    known_pos_range2& localEvidenceRange)
{
    const SimpleAlignment localAlign(localRead);

    // 1) process any large indels in the localRead:
    getSVCandidatesFromRead(opt, localAlign, localRead.target_id(), candidates);

    // 2) process anomalous read pair relationships:
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

        localEvidenceRange.set_range(startRefPos,endRefPos);
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

    // check if read pair separation is non-anomalous after accounting for read alignments:
    if ((localRead.target_id() == localRead.mate_target_id()) &&
        (localRead.is_fwd_strand() != localRead.is_mate_fwd_strand()))
    {
        // get length of fragment after accounting for any variants described directly in either read alignment:
        const pos_t cigarAdjustedFragmentSize(totalNoninsertSize + (insertRange.end_pos() - insertRange.begin_pos()));
        if ((cigarAdjustedFragmentSize <= rstats.properPair.max) &&
            (cigarAdjustedFragmentSize >= rstats.properPair.min)) return;
    }

    candidates.push_back(sv);
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
    const SVLocusScanner::CachedReadGroupStats& rstats,
    const bam_record& bamRead,
    std::vector<SVLocus>& loci)
{
    using namespace illumina::common;

    loci.clear();
    std::vector<SVCandidate> candidates;
    known_pos_range2 localEvidenceRange;

    getReadBreakendsImpl(opt, rstats, bamRead, NULL, candidates, localEvidenceRange);

    BOOST_FOREACH(const SVCandidate& cand, candidates)
    {
        const SVBreakend& localBreakend(cand.bp1);
        const SVBreakend& remoteBreakend(cand.bp2);
        if ((0==localBreakend.interval.range.size()) ||
            (0==remoteBreakend.interval.range.size()))
        {
            std::ostringstream oss;
            oss << "Empty Breakend proposed from bam record.\n"
                << "\tlocal_breakend: " << localBreakend << "\n"
                << "\tremote_breakend: " << remoteBreakend << "\n"
                << "\tbam_record: " << bamRead << "\n";
            BOOST_THROW_EXCEPTION(LogicException(oss.str()));
        }

        SVLocus locus;
        // set local breakend estimate:
        NodeIndexType localBreakendNode(0);
        {
            localBreakendNode = locus.addNode(localBreakend.interval);
            locus.setNodeEvidence(localBreakendNode,localEvidenceRange);
        }

        // set remote breakend estimate:
        {
            NodeIndexType remoteBreakendNode;
            if ((remoteBreakend.readCount != 0) ||
                (remoteBreakend.splitCount != 0))
            {
                remoteBreakendNode = locus.addNode(remoteBreakend.interval);
                locus.linkNodes(localBreakendNode,remoteBreakendNode,1,1);
            }
            else
            {
                remoteBreakendNode = locus.addRemoteNode(remoteBreakend.interval);
                locus.linkNodes(localBreakendNode,remoteBreakendNode);
            }
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
    _opt(opt)
{
    // pull in insert stats:
    _rss.read(statsFilename.c_str());

    // cache the insert stats we'll be looking up most often:
    BOOST_FOREACH(const std::string& file, alignmentFilename)
    {
        const boost::optional<unsigned> index(_rss.getGroupIndex(file));
        assert(index);
        const ReadGroupStats rgs(_rss.getStats(*index));

        _stats.resize(_stats.size()+1);
        CachedReadGroupStats& stat(_stats.back());
        {
            Range& breakend(stat.breakendRegion);
            breakend.min=rgs.fragSize.quantile(_opt.breakendEdgeTrimProb);
            breakend.max=rgs.fragSize.quantile((1-_opt.breakendEdgeTrimProb));

            if (breakend.min<0.) breakend.min = 0;
            assert(breakend.max>0.);
        }
        {
            Range& ppair(stat.properPair);
            ppair.min=rgs.fragSize.quantile(_opt.properPairTrimProb);
            ppair.max=rgs.fragSize.quantile((1-_opt.properPairTrimProb));

            if (ppair.min<0.) ppair.min = 0;

            assert(ppair.max>0.);
        }
        {
            Range& evid(stat.evidencePair);
            evid.min=rgs.fragSize.quantile(_opt.evidenceTrimProb);
            evid.max=rgs.fragSize.quantile((1-_opt.evidenceTrimProb));

            if (evid.min<0.) evid.min = 0;

            assert(evid.max>0.);
        }
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
    if (bamRead.is_unmapped() || bamRead.is_mate_unmapped()) return false;
    if (bamRead.target_id() != bamRead.mate_target_id()) return false;

    const Range& ppr(_stats[defaultReadGroupIndex].properPair);
    const int32_t fragmentSize(std::abs(bamRead.template_size()));
    if (fragmentSize > ppr.max || fragmentSize < ppr.min) return false;

    if     (bamRead.pos() < bamRead.mate_pos())
    {
        if (! bamRead.is_fwd_strand()) return false;
        if (  bamRead.is_mate_fwd_strand()) return false;
    }
    else if (bamRead.pos() > bamRead.mate_pos())
    {
        if (  bamRead.is_fwd_strand()) return false;
        if (! bamRead.is_mate_fwd_strand()) return false;
    }
    else
    {
        if (bamRead.is_fwd_strand() == bamRead.is_mate_fwd_strand()) return false;
    }

    return true;
}



void
SVLocusScanner::
getSVLoci(
    const bam_record& bamRead,
    const unsigned defaultReadGroupIndex,
    std::vector<SVLocus>& loci) const
{
    loci.clear();

    if (! bamRead.is_chimeric())
    {
        if (std::abs(bamRead.template_size())<2000) return;
    }

    const CachedReadGroupStats& rstats(_stats[defaultReadGroupIndex]);
    getSVLociImpl(_opt, rstats, bamRead, loci);
}



void
SVLocusScanner::
getBreakendPair(
    const bam_record& localRead,
    const bam_record* remoteReadPtr,
    const unsigned defaultReadGroupIndex,
    std::vector<SVCandidate>& candidates) const
{
    const CachedReadGroupStats& rstats(_stats[defaultReadGroupIndex]);

    // throw evidence range away in this case
    known_pos_range2 evidenceRange;
    getReadBreakendsImpl(_opt, rstats, localRead, remoteReadPtr, candidates, evidenceRange);
}
