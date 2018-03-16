//
// Manta - Structural Variant and Indel Caller
// Copyright (c) 2013-2018 Illumina, Inc.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//

/// \file
/// \author Chris Saunders
///

#include "SVLocusSetFinder.hh"

#include "blt_util/log.hh"
#include "common/Exceptions.hh"
#include "htsapi/align_path_bam_util.hh"
#include "manta/ChromDepthFilterUtil.hh"
#include "manta/ReadFilter.hh"

#include <iostream>
#include <sstream>


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
    // the number of positions below the HEAD position where the
    // depth buffer can be safely erased
    static const unsigned clearDepthBorderSize(10);

    stage_data sd;
    sd.add_stage(HEAD);
    sd.add_stage(DENOISE, HEAD, denoiseBorderSize);
    sd.add_stage(CLEAR_DEPTH, HEAD, clearDepthBorderSize);

    return sd;
}
}



/// \brief Compute the genomic region in which graph denoising is allowed.
///
/// The allowed denoising region is the same as the current scan region, except that the region is shortened
/// by a fixed-size buffer wherever an adjacent genome segment is potentially being analyzed by another
/// process in parallel. This protected zone near adjacent segment boundaries helps to ensure that SV
/// breakend regions spanning the adjacent process boundary are not pessimistically marked as noise before
/// the evidence from both segments has been merged.
///
/// \param scanRegion The region of the genome scanned by this process for SV locus evidence.
/// \param bamHeader Bam header information. Used to extract chromosome length here.
/// \param denoiseRegionProtectedBorderSize Length of the protected region where denoising is skipped at adjacent
///                                         process boundaries.
/// \return Region where denoising is allowed in this process
static
GenomeInterval
computeDenoiseRegion(
    const GenomeInterval& scanRegion,
    const bam_header_info& bamHeader,
    const int denoiseRegionProtectedBorderSize)
{
    GenomeInterval denoiseRegion=scanRegion;

    known_pos_range2& range(denoiseRegion.range);
    if (range.begin_pos() > 0)
    {
        range.set_begin_pos(range.begin_pos()+denoiseRegionProtectedBorderSize);
    }

    bool isEndBorder(true);
    if (static_cast<int32_t>(bamHeader.chrom_data.size()) > denoiseRegion.tid)
    {
        const pos_t chromEndPos(bamHeader.chrom_data[denoiseRegion.tid].length);
        isEndBorder=(range.end_pos() < chromEndPos);
    }

    if (isEndBorder)
    {
        range.set_end_pos(range.end_pos()-denoiseRegionProtectedBorderSize);
    }

#ifdef DEBUG_SFINDER
    log_os << __FUNCTION__ << ": " << denoiseRegion << "\n";
#endif

    return  denoiseRegion;
}


/// The compression level used by this object's private depth-per-position buffer.
static const unsigned depthBufferCompression = 16;



SVLocusSetFinder::
SVLocusSetFinder(
    const ESLOptions& opt,
    const GenomeInterval& scanRegion,
    const std::shared_ptr<reference_contig_segment> refSeqPtr,
    std::shared_ptr<SVLocusSet> svLociPtr) :
    _isAlignmentTumor(opt.alignFileOpt.isAlignmentTumor),
    _scanRegion(scanRegion),
    _refSeqPtr(refSeqPtr),
    _svLociPtr(svLociPtr),
    _denoiseRegion(computeDenoiseRegion(scanRegion, _bamHeader(), REGION_DENOISE_BORDER)),
    _stageManager(
        STAGE::getStageData(REGION_DENOISE_BORDER),
        pos_range(
            scanRegion.range.begin_pos(),
            scanRegion.range.end_pos()),
        *this),
    _positionReadDepthEstimate(depthBufferCompression),
    _isInDenoiseRegion(false),
    _denoiseStartPos(0),
    _readScanner(opt.scanOpt,opt.statsFilename,opt.alignFileOpt.alignmentFilenames, opt.isRNA),
    _isMaxDepthFilter(false),
    _maxDepth(0)
{
    assert(_svLociPtr);

    //
    // initialize max depth filtration values:
    //
    const ChromDepthFilterUtil dFilter(opt.chromDepthFilename, opt.scanOpt.maxDepthFactor, _bamHeader());
    _isMaxDepthFilter=dFilter.isMaxDepthFilter();
    if (_isMaxDepthFilter)
    {
        _maxDepth=dFilter.maxDepth(scanRegion.tid);
    }

    // If SV locus graph is empty, initialize various metadata, otherwise verify that meta-data match existing metadata
    //
    const unsigned sampleCount(opt.alignFileOpt.alignmentFilenames.size());
    if (getLocusSet().getCounts().size() == 0)
    {
        _getLocusSet().getCounts().setSampleCount(sampleCount);
        for (unsigned sampleIndex(0); sampleIndex < sampleCount; ++sampleIndex)
        {
            _getLocusSet().getCounts().getSampleCounts(sampleIndex).sampleSource =
                opt.alignFileOpt.alignmentFilenames[sampleIndex];
        }
    }
    else
    {
        assert(getLocusSet().getCounts().size() == sampleCount);
        for (unsigned sampleIndex(0); sampleIndex < sampleCount; ++sampleIndex)
        {
            assert(getLocusSet().getCounts().getSampleCounts(sampleIndex).sampleSource ==
                       opt.alignFileOpt.alignmentFilenames[sampleIndex]);
        }
    }
}



void
SVLocusSetFinder::
process_pos(const int stage_no,
            const pos_t pos)
{
#ifdef DEBUG_SFINDER
    log_os << __FUNCTION__ << ": stage_no: " << stage_no << " pos: " << pos << "\n";
#endif

    if     (stage_no == STAGE::HEAD)
    {
        // pass
    }
    else if (stage_no == STAGE::DENOISE)
    {
        // denoise the SV locus graph in regions of at least this size
        //
        // this parameter batches the denoising process to balance runtime vs.
        // maintaining the SV locus graph as a reasonably compact data structure
        //
        static const pos_t minDenoiseRegionSize(1000);

        if (_denoiseRegion.range.is_pos_intersect(pos))
        {

#ifdef DEBUG_SFINDER
            log_os << __FUNCTION__ << ": pos intersect. pos: " << pos << " dnRegion: " << _denoiseRegion << " is in region: " << _isInDenoiseRegion << "\n";
#endif

            if (! _isInDenoiseRegion)
            {
                _denoiseStartPos=_denoiseRegion.range.begin_pos();
                _isInDenoiseRegion=true;
            }

            if ( (1 + pos-_denoiseStartPos) >= minDenoiseRegionSize)
            {
                _getLocusSet().cleanRegion(GenomeInterval(_denoiseRegion.tid, _denoiseStartPos, (pos+1)));
                _denoiseStartPos = (pos+1);
            }
        }
        else
        {

#ifdef DEBUG_SFINDER
            log_os << __FUNCTION__ << ": no pos intersect. pos: " << pos << " dnRegion: " << _denoiseRegion << " is in region: " << _isInDenoiseRegion << "\n";
#endif

            if (_isInDenoiseRegion)
            {
                if ( (_denoiseRegion.range.end_pos()-_denoiseStartPos) > 0)
                {
                    _getLocusSet().cleanRegion(GenomeInterval(_denoiseRegion.tid, _denoiseStartPos, _denoiseRegion.range.end_pos()));
                    _denoiseStartPos = _denoiseRegion.range.end_pos();
                }
                _isInDenoiseRegion=false;
            }
        }
    }
    else if (stage_no == STAGE::CLEAR_DEPTH)
    {
        _positionReadDepthEstimate.clear_pos(pos);
    }
    else
    {
        assert(false && "Unexpected stage id");
    }
}



void
SVLocusSetFinder::
addToDepthBuffer(
    const bam_record& bamRead)
{
    if (! _isMaxDepthFilter) return;

    const pos_t refPos(bamRead.pos()-1);

    // Estimated read depth uses a very simple approximation that each input read aligns without any indels.
    // This is done to reduce the depth estimation runtime overhead.
    const pos_t readSize(bamRead.read_size());
    _positionReadDepthEstimate.inc(refPos,readSize);
}



void
SVLocusSetFinder::
update(
    const stream_state_reporter& streamErrorReporter,
    const bam_record& bamRead,
    const unsigned defaultReadGroupIndex)
{
    // This is the primary read filtration step for the purpose of graph building
    //
    // Although unmapped reads themselves are filtered out, reads with unmapped mates are still
    // accepted, because these contribute signal for assembly (indel) regions.
    if (isReadUnmappedOrFilteredCore(bamRead)) return;

    const bool isTumor(_isAlignmentTumor[defaultReadGroupIndex]);
    if (! isTumor) addToDepthBuffer(bamRead);

    // Filter out reads from high-depth chromosome regions
    if (_isMaxDepthFilter)
    {
        if (_positionReadDepthEstimate.val(bamRead.pos()-1) > _maxDepth) return;
    }

    // Verify the input reads conform to BAM input restrictions before testing if they could be input evidence
    // The reason this is done here is that the restriction against the '=' sign in the sequence could otherwise
    // silently alter downstream evidence checks and assembly steps:
    {
        const bam_seq readSeq(bamRead.get_bam_read());
        const unsigned readSize(readSeq.size());
        for (unsigned baseIndex(0); baseIndex < readSize; baseIndex++)
        {
            if (readSeq.get_code(baseIndex) == BAM_BASE::REF)
            {
                std::ostringstream oss;
                oss << "Unsupported use of the '=' symbol in the BAM/CRAM SEQ field from read:\n";
                streamErrorReporter.report_state(oss);
                BOOST_THROW_EXCEPTION(illumina::common::GeneralException(oss.str()));
            }
        }
    }

    // counts/incounts are part of the statistics tracking framework used for
    // methods diagnostics. They do not impact the graph build.
    SampleCounts& counts(_getLocusSet().getCounts().getSampleCounts(defaultReadGroupIndex));
    SampleReadInputCounts& incounts(counts.input);

    // Filter out reads below the minimum MAPQ threshold
    if (bamRead.map_qual() < _readScanner.getMinMapQ())
    {
        incounts.minMapq++;
        return;
    }


    // Filter out reads which are not found to be indicative of an SV
    //
    // For certain types of SV evidence, these tests are designed to be a faster approximation of the
    // full test which will be repeated when the read(-pairs) are converted to SV locus graph elements,
    // therefor a small number of reads which make it through this filter will be effectively filtered
    // later as non-evidence -- they will not be filtered per se, but rather converted into zero SVLocus
    // objects.
    //
    if (! _readScanner.isSVEvidence(bamRead,defaultReadGroupIndex,_refSeq(),&(incounts.evidenceCount))) return;

#ifdef DEBUG_SFINDER
    log_os << __FUNCTION__ << ": Accepted read. isNonCompressedAnomalous "  << isNonCompressedAnomalous << " is Local assm evidence: " << isLocalAssemblyEvidence << " read: " << bamRead << "\n";
#endif

    // check that this read starts in our scan region:
    if (! _scanRegion.range.is_pos_intersect(bamRead.pos()-1)) return;

    // QC check of read length
    SVLocusScanner::checkReadSize(streamErrorReporter, bamRead);

    // update the stage manager to move the head pointer forward to the current read's position
    //
    // This may trigger a denoising step or clear buffered information for all positions at a
    // given offset below the new head position
    //
    _stageManager.handle_new_pos_value(bamRead.pos()-1);

    // convert the given read into zero to many SVLocus objects
    //
    // In almost all cases, each read should be converted into one SVLocus object, and that
    // SVLocus will consist of either one or two SVLocusNodes.
    //
    std::vector<SVLocus> loci;
    SampleEvidenceCounts& eCounts(counts.evidence);
    _readScanner.getSVLoci(bamRead, defaultReadGroupIndex, _bamHeader(),
                           _refSeq(), loci, eCounts);

    // merge each non-empty SV locus into this genome segment graph:
    for (const SVLocus& locus : loci)
    {
        if (locus.empty()) continue;
        _getLocusSet().merge(locus);
    }
}
