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
    }
}



void
SVLocusScanner::
getChimericSVLocusImpl(
    const CachedReadGroupStats& rstats,
    const bam_record& read,
    SVLocus& locus)
{
    ALIGNPATH::path_t apath;
    bam_cigar_to_apath(read.raw_cigar(),read.n_cigar(),apath);

    const unsigned readSize(apath_read_length(apath));

    unsigned thisReadNoninsertSize(0);
    if (read.is_fwd_strand())
    {
        thisReadNoninsertSize=(readSize-apath_read_trail_size(apath));
    }
    else
    {
        thisReadNoninsertSize=(readSize-apath_read_lead_size(apath));
    }

    // estimate mate read size to be same as this, and no clipping on mate read:
    const unsigned mateReadNoninsertSize(readSize);

    const unsigned totalNoninsertSize(thisReadNoninsertSize+mateReadNoninsertSize);

    // get local breakend estimate:
    NodeIndexType localBreakendNode(0);
    {
        GenomeInterval localBreakend(read.target_id());

        const pos_t startRefPos(read.pos()-1);
        const pos_t endRefPos(startRefPos+apath_ref_length(apath));
        // expected breakpoint range is from the end of the read alignment to the (probabilistic) end of the fragment:
        if (read.is_fwd_strand())
        {
            localBreakend.range.set_begin_pos(endRefPos);
            localBreakend.range.set_end_pos(endRefPos + static_cast<pos_t>(rstats.max-(totalNoninsertSize)));
        }
        else
        {
            localBreakend.range.set_end_pos(startRefPos);
            localBreakend.range.set_begin_pos(startRefPos - static_cast<pos_t>(rstats.max-(totalNoninsertSize)));
        }

        const known_pos_range2 evidenceRange(startRefPos,endRefPos);
        localBreakendNode = locus.addNode(localBreakend);
        locus.setNodeEvidence(localBreakendNode,evidenceRange);
    }

    // get remote breakend estimate:
    {
        GenomeInterval remoteBreakend(read.mate_target_id());

        const pos_t startRefPos(read.mate_pos()-1);
        const pos_t endRefPos(startRefPos+readSize);
        if (read.is_mate_fwd_strand())
        {
            remoteBreakend.range.set_begin_pos(endRefPos);
            remoteBreakend.range.set_end_pos(endRefPos + static_cast<pos_t>(rstats.max-(totalNoninsertSize)));
        }
        else
        {
            remoteBreakend.range.set_end_pos(startRefPos);
            remoteBreakend.range.set_begin_pos(startRefPos - static_cast<pos_t>(rstats.max-(totalNoninsertSize)));
        }

        const NodeIndexType remoteBreakendNode(locus.addRemoteNode(remoteBreakend));
        locus.linkNodes(localBreakendNode,remoteBreakendNode);
    }

}



bool
SVLocusScanner::
isReadFiltered(const bam_record& read) const
{
    if (read.is_filter()) return true;
    if (read.is_dup()) return true;
    if (read.is_secondary()) return true;
    if (read.is_proper_pair()) return true;
    if (read.map_qual() < _opt.minMapq) return true;
    return false;
}



void
SVLocusScanner::
getChimericSVLocus(const bam_record& read,
                   const unsigned defaultReadGroupIndex,
                   SVLocus& locus) const
{
    locus.clear();

    if (read.is_chimeric())
    {
        const CachedReadGroupStats& rstats(_stats[defaultReadGroupIndex]);
        getChimericSVLocusImpl(rstats,read,locus);
    }
}
