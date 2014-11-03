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
///

#include "SVLocusSetFinder.hh"

#include "blt_util/align_path_bam_util.hh"
#include "blt_util/log.hh"
#include "manta/ChromDepthFilterUtil.hh"

#include <iostream>



namespace STAGE
{
enum index_t
{
    HEAD,
    DENOISE,
    CLEAR_DEPTH
};


static
stage_data
getStageData(
    const unsigned denoiseBorderSize)
{
    static const unsigned clearDepthBorderSize(10);

    stage_data sd;
    sd.add_stage(HEAD);
    sd.add_stage(DENOISE, HEAD, denoiseBorderSize);
    sd.add_stage(CLEAR_DEPTH, HEAD, clearDepthBorderSize);

    return sd;
}
}



SVLocusSetFinder::
SVLocusSetFinder(
    const ESLOptions& opt,
    const GenomeInterval& scanRegion,
    const bam_header_info& bamHeader) :
    _isAlignmentTumor(opt.alignFileOpt.isAlignmentTumor),
    _scanRegion(scanRegion),
    _stageman(
        STAGE::getStageData(REGION_DENOISE_BORDER),
        pos_range(
            scanRegion.range.begin_pos(),
            scanRegion.range.end_pos()),
        *this),
    _svLoci(opt.graphOpt),
    _isScanStarted(false),
    _isInDenoiseRegion(false),
    _denoisePos(0),
    _readScanner(opt.scanOpt,opt.statsFilename,opt.alignFileOpt.alignmentFilename),
    _isMaxDepth(false),
    _maxDepth(0)
{
    const ChromDepthFilterUtil dFilter(opt.chromDepthFilename, opt.scanOpt.maxDepthFactor, bamHeader);
    _isMaxDepth=dFilter.isMaxDepthFilter();
    if (_isMaxDepth)
    {
        _maxDepth=dFilter.maxDepth(scanRegion.tid);
    }

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
            log_os << "SFinder::process_pos pos intersect. pos: " << pos << " dnRegion: " << _denoiseRegion << " is in region: " << _isInDenoiseRegion << "\n";
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
            log_os << "SFinder::process_pos no pos intersect. pos: " << pos << " dnRegion: " << _denoiseRegion << " is in region: " << _isInDenoiseRegion << "\n";
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
    else if (stage_no == STAGE::CLEAR_DEPTH)
    {
        _depth.clear_pos(pos);
    }
    else
    {
        assert(false && "Unexpected stage id");
    }
}



void
SVLocusSetFinder::
addToDepthBuffer(
    const unsigned defaultReadGroupIndex,
    const bam_record& bamRead)
{
    if (! _isMaxDepth) return;

    // estimate depth from normal sample only:
    if (_isAlignmentTumor[defaultReadGroupIndex]) return;

    // depth estimation relies on a simple filtration criteria to stay in sync with the chromosome mean
    // depth estimates:
    if (bamRead.is_unmapped()) return;

    const pos_t refPos(bamRead.pos()-1);

    /// stick to a simple approximation -- ignore CIGAR string and just look at the read length:
    const pos_t readSize(bamRead.read_size());
    for (pos_t readIndex(0); readIndex<readSize; ++readIndex)
    {
        _depth.inc(refPos+readIndex);
    }
}



void
SVLocusSetFinder::
update(
    const bam_record& bamRead,
    const unsigned defaultReadGroupIndex,
    const bam_header_info& bamHeader,
    const reference_contig_segment& refSeq,
    TruthTracker& truthTracker)
{
    _isScanStarted=true;

    const bool isTumor(_isAlignmentTumor[defaultReadGroupIndex]);
    if (! isTumor)
    {
        // depth estimation relies on a simple filtration criteria to stay in sync with the chromosome mean
        // depth estimates:
        if (! bamRead.is_unmapped())
        {
            addToDepthBuffer(defaultReadGroupIndex, bamRead);
        }
    }

    if (_readScanner.isReadFiltered(bamRead)) return;

    if (_isMaxDepth)
    {
        if (_depth.val(bamRead.pos()-1) > _maxDepth) return;
    }

    // exclude innie read pairs which are anomalously short:
    const bool isNonCompressedAnomalous(_readScanner.isNonCompressedAnomalous(bamRead,defaultReadGroupIndex));

    if (isNonCompressedAnomalous) _svLoci.getReadCounts(isTumor).anom++;
    else                          _svLoci.getReadCounts(isTumor).nonAnom++;

    bool isLocalAssemblyEvidence(false);
    if (! isNonCompressedAnomalous)
    {
        isLocalAssemblyEvidence = _readScanner.isLocalAssemblyEvidence(bamRead, refSeq);
    }

    const bool isRejectRead(! ( isNonCompressedAnomalous || isLocalAssemblyEvidence));

    if (isRejectRead)
    {
        return; // this read isn't interesting (enough) wrt SV discovery
    }

#ifdef DEBUG_SFINDER
    log_os << "SFinder: Accepted read. isNonCompressedAnomalous "  << isNonCompressedAnomalous << " is Local assm evidence: " << isLocalAssemblyEvidence << " read: " << bamRead << "\n";
#endif

    // check that this read starts in our scan region:
    if (! _scanRegion.range.is_pos_intersect(bamRead.pos()-1)) return;

    _stageman.handle_new_pos_value(bamRead.pos()-1);

    std::vector<SVLocus> loci;

    _readScanner.getSVLoci(bamRead, defaultReadGroupIndex, bamHeader,
                           refSeq, loci, truthTracker);

    for (const SVLocus& locus : loci)
    {
        if (locus.empty()) continue;
        _svLoci.merge(locus);
    }
}
