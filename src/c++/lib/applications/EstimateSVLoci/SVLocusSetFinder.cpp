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

#include "SVLocusSetFinder.hh"
#include "blt_util/align_path_bam_util.hh"
#include "blt_util/log.hh"

#include "boost/foreach.hpp"

#include <iostream>



namespace STAGE {
    enum index_t {
        HEAD,
        DENOISE
    };


    static
    stage_data
    getStageData(const unsigned denoiseBorderSize)
    {
        stage_data sd;
        sd.add_stage(HEAD);
        sd.add_stage(DENOISE, HEAD, denoiseBorderSize);

        return sd;
    }
}



SVLocusSetFinder::
SVLocusSetFinder(
        const ESLOptions& opt,
        const GenomeInterval& scanRegion) :
        _opt(opt),
        _scanRegion(scanRegion),
        _stageman(
                STAGE::getStageData(REGION_DENOISE_BORDER),
                pos_range(
                        scanRegion.range.begin_pos(),
                        scanRegion.range.end_pos()),
                *this),
        _isScanStarted(false),
        _isInDenoiseRegion(false),
        _denoisePos(0)
{
    // pull in insert stats:
    _rss.read(opt.statsFilename.c_str());

    updateDenoiseRegion();

    // cache the insert stats we'll be looking up most often:
    BOOST_FOREACH(const std::string& file, opt.alignmentFilename)
    {
        const boost::optional<unsigned> index(_rss.getGroupIndex(file));
        assert(index);
        const ReadGroupStats rgs(_rss.getStats(*index));

        _stats.resize(_stats.size()+1);
        _stats.back().min=rgs.fragSize.quantile(opt.breakendEdgeTrimProb);
        _stats.back().max=rgs.fragSize.quantile((1-opt.breakendEdgeTrimProb));
    }
}



void
SVLocusSetFinder::
updateDenoiseRegion()
{
    _denoiseRegion=_scanRegion;

    known_pos_range2& range(_denoiseRegion.range);
    if (range.begin_pos() > 0)
    {
        range.set_begin_pos(range.begin_pos()+REGION_DENOISE_BORDER);
    }

    bool isEndBorder(true);
    if(static_cast<int32_t>(_svLoci.header.chrom_data.size()) > _denoiseRegion.tid)
    {
        const pos_t chromEndPos(_svLoci.header.chrom_data[_denoiseRegion.tid].length);
        isEndBorder=(range.end_pos() < chromEndPos);
    }

    if (isEndBorder)
    {
        range.set_end_pos(range.end_pos()-REGION_DENOISE_BORDER);
    }

#ifdef DEBUG_SFINDER
    log_os << "SFinder::updateDenoiseRegion " << _denoiseRegion << "\n";
#endif

}



void
SVLocusSetFinder::
process_pos(const int stage_no,
            const pos_t pos)
{

#ifdef DEBUG_SFINDER
    log_os << "SFinder::process_pos stage_no: " << stage_no << " pos: " << pos << "\n";
#endif

    if     (stage_no == STAGE::HEAD)
    {
        // pass
    }
    else if(stage_no == STAGE::DENOISE)
    {
        static const pos_t denoiseMinChunk(1000);

        if(_denoiseRegion.range.is_pos_intersect(pos))
        {

#ifdef DEBUG_SFINDER
    log_os << "SFinder::process_pos pos intersect. is in region: " << _isInDenoiseRegion << "\n";
#endif

            if(! _isInDenoiseRegion)
            {
                _denoisePos=_denoiseRegion.range.begin_pos();
                _isInDenoiseRegion=true;
            }

            if( (1 + pos-_denoisePos) >= denoiseMinChunk)
            {
                _svLoci.cleanRegion(GenomeInterval(_denoiseRegion.tid, _denoisePos, (pos+1)));
                _denoisePos = (pos+1);
            }
        }
        else
        {

#ifdef DEBUG_SFINDER
    log_os << "SFinder::process_pos pos intersect. is in region: " << _isInDenoiseRegion << "\n";
#endif

            if(_isInDenoiseRegion)
            {
                if( (_denoiseRegion.range.end_pos()-_denoisePos) > 0)
                {
                    _svLoci.cleanRegion(GenomeInterval(_denoiseRegion.tid, _denoisePos, _denoiseRegion.range.end_pos()));
                    _denoisePos = _denoiseRegion.range.end_pos();
                }
                _isInDenoiseRegion=false;
            }
        }
    }
    else
    {
        assert(0);
    }
}



void
SVLocusSetFinder::
getChimericSVLocus(
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
SVLocusSetFinder::
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
SVLocusSetFinder::
update(const bam_record& read,
       const unsigned defaultReadGroupIndex)
{
    _isScanStarted=true;

    // shortcut to speed things up:
    if (isReadFiltered(read)) return;

    // check that this read starts in our scan region:
    if(! _scanRegion.range.is_pos_intersect(read.pos()-1)) return;

    _stageman.handle_new_pos_value(read.pos()-1);

    const CachedReadGroupStats& rstats(_stats[defaultReadGroupIndex]);

    SVLocus locus;

    // start out looking for chimeric reads only:
    if (read.is_chimeric())
    {
        getChimericSVLocus(rstats,read,locus);
    }

    if (! locus.empty())
    {
        _svLoci.merge(locus);
    }
}
