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



namespace STAGE
{
enum index_t
{
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
    _scanRegion(scanRegion),
    _stageman(
        STAGE::getStageData(REGION_DENOISE_BORDER),
        pos_range(
            scanRegion.range.begin_pos(),
            scanRegion.range.end_pos()),
        *this),
    _svLoci(opt.minMergeEdgeCount),
    _isScanStarted(false),
    _isInDenoiseRegion(false),
    _denoisePos(0),
    _readScanner(opt.scanOpt,opt.statsFilename,opt.alignmentFilename),
    _anomCount(0),
    _nonAnomCount(0)
{
    updateDenoiseRegion();
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
    if (static_cast<int32_t>(_svLoci.header.chrom_data.size()) > _denoiseRegion.tid)
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
    else if (stage_no == STAGE::DENOISE)
    {
        static const pos_t denoiseMinChunk(1000);

        if (_denoiseRegion.range.is_pos_intersect(pos))
        {

#ifdef DEBUG_SFINDER
            log_os << "SFinder::process_pos pos intersect. is in region: " << _isInDenoiseRegion << "\n";
#endif

            if (! _isInDenoiseRegion)
            {
                _denoisePos=_denoiseRegion.range.begin_pos();
                _isInDenoiseRegion=true;
            }

            if ( (1 + pos-_denoisePos) >= denoiseMinChunk)
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

            if (_isInDenoiseRegion)
            {
                if ( (_denoiseRegion.range.end_pos()-_denoisePos) > 0)
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
update(const bam_record& bamRead,
       const unsigned defaultReadGroupIndex)
{
    _isScanStarted=true;

    // shortcut to speed things up:
    if (_readScanner.isReadFiltered(bamRead)) return;
    if (bamRead.is_proper_pair())
    {
        _nonAnomCount++;
        return;
    }
    _anomCount++;

    // check that this read starts in our scan region:
    if (! _scanRegion.range.is_pos_intersect(bamRead.pos()-1)) return;

    _stageman.handle_new_pos_value(bamRead.pos()-1);

    SVLocus locus;

    //_readScanner.getSVLocus(bamRead, defaultReadGroupIndex, locus);
    _readScanner.getSVLocus(bamRead, defaultReadGroupIndex, locus);

    if (! locus.empty())
    {
        _svLoci.merge(locus);
    }
}
