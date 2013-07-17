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
#include "blt_util/log.hh"
#include "common/Exceptions.hh"

#include "boost/foreach.hpp"



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
        _stats.back().min=rgs.fragSize.quantile(_opt.breakendEdgeTrimProb);
        _stats.back().max=rgs.fragSize.quantile((1-_opt.breakendEdgeTrimProb));

        if(_stats.back().min<0.) _stats.back().min = 0;

        assert(_stats.back().max>0.);
    }
}



void
SVLocusScanner::
getReadBreakendsImpl(
    const CachedReadGroupStats& rstats,
    const bam_record& localRead,
    const bam_record* remoteReadPtr,
    SVBreakend& localBreakend,
    SVBreakend& remoteBreakend,
    known_pos_range2& evidenceRange)
{
    static const pos_t minPairBreakendSize(40);

    ALIGNPATH::path_t apath;
    bam_cigar_to_apath(localRead.raw_cigar(),localRead.n_cigar(),apath);

    const unsigned readSize(apath_read_length(apath));
    const unsigned localRefLength(apath_ref_length(apath));

    unsigned thisReadNoninsertSize(0);
    if (localRead.is_fwd_strand())
    {
        thisReadNoninsertSize=(readSize-apath_read_trail_size(apath));
    }
    else
    {
        thisReadNoninsertSize=(readSize-apath_read_lead_size(apath));
    }

    localBreakend.readCount = 1;

    // if remoteRead is not available, estimate mate localRead size to be same as local,
    // and assume no clipping on mate localRead:
    unsigned remoteReadNoninsertSize(readSize);
    unsigned remoteRefLength(localRefLength);

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

    const pos_t totalNoninsertSize(thisReadNoninsertSize+remoteReadNoninsertSize);
    const pos_t breakendSize(std::max(minPairBreakendSize,static_cast<pos_t>(rstats.max-totalNoninsertSize)));

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
        }
        else
        {
            localBreakend.state = SVBreakendState::LEFT_OPEN;
            localBreakend.interval.range.set_end_pos(startRefPos);
            localBreakend.interval.range.set_begin_pos(startRefPos - breakendSize);
        }

        evidenceRange.set_range(startRefPos,endRefPos);
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
        }
        else
        {
            remoteBreakend.state = SVBreakendState::LEFT_OPEN;
            remoteBreakend.interval.range.set_end_pos(startRefPos);
            remoteBreakend.interval.range.set_begin_pos(startRefPos - breakendSize);
        }
    }
}



void
SVLocusScanner::
getSVLocusImpl(
    const CachedReadGroupStats& rstats,
    const bam_record& read,
    SVLocus& locus)
{
    using namespace illumina::common;

    SVBreakend localBreakend;
    SVBreakend remoteBreakend;
    known_pos_range2 evidenceRange;
    getReadBreakendsImpl(rstats, read, NULL, localBreakend, remoteBreakend, evidenceRange);

    if((0==localBreakend.interval.range.size()) ||
       (0==remoteBreakend.interval.range.size()))
    {
        std::ostringstream oss;
        oss << "Empty Breakend proposed from bam record.\n"
            << "\tlocal_breakend: " << localBreakend << "\n"
            << "\tremote_breakend: " << remoteBreakend << "\n"
            << "\tbam_record: " << read << "\n";
        BOOST_THROW_EXCEPTION(LogicException(oss.str()));
    }

    // set local breakend estimate:
    NodeIndexType localBreakendNode(0);
    {
        localBreakendNode = locus.addNode(localBreakend.interval);
        locus.setNodeEvidence(localBreakendNode,evidenceRange);
    }

    // set remote breakend estimate:
    {
        const NodeIndexType remoteBreakendNode(locus.addRemoteNode(remoteBreakend.interval));
        locus.linkNodes(localBreakendNode,remoteBreakendNode);
        locus.mergeSelfOverlap();
    }

}



bool
SVLocusScanner::
isReadFiltered(const bam_record& read) const
{
    if (read.is_filter()) return true;
    if (read.is_dup()) return true;
    if (read.is_secondary()) return true;
    if (read.map_qual() < _opt.minMapq) return true;
    return false;
}



void
SVLocusScanner::
getChimericSVLocus(
    const bam_record& read,
    const unsigned defaultReadGroupIndex,
    SVLocus& locus) const
{
    locus.clear();

    if (read.is_chimeric())
    {
        const CachedReadGroupStats& rstats(_stats[defaultReadGroupIndex]);
        getSVLocusImpl(rstats,read,locus);
    }
}



void
SVLocusScanner::
getSVLocus(
    const bam_record& read,
    const unsigned defaultReadGroupIndex,
    SVLocus& locus) const
{
    locus.clear();

    if(! read.is_chimeric())
    {
        if(read.template_size()<1000) return;
    }

    const CachedReadGroupStats& rstats(_stats[defaultReadGroupIndex]);
    getSVLocusImpl(rstats,read,locus);
}



void
SVLocusScanner::
getBreakendPair(
    const bam_record& localRead,
    const bam_record* remoteReadPtr,
    const unsigned defaultReadGroupIndex,
    SVBreakend& localBreakend,
    SVBreakend& remoteBreakend) const
{
    const CachedReadGroupStats& rstats(_stats[defaultReadGroupIndex]);
    known_pos_range2 evidenceRange;
    getReadBreakendsImpl(rstats, localRead, remoteReadPtr, localBreakend, remoteBreakend, evidenceRange);
}
