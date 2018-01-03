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

#pragma once

#include "ESLOptions.hh"

#include "blt_util/depth_buffer.hh"
#include "blt_util/pos_processor_base.hh"
#include "blt_util/stage_manager.hh"
#include "htsapi/bam_record.hh"
#include "htsapi/bam_streamer.hh"
#include "manta/SVLocusScanner.hh"
#include "svgraph/SVLocusSet.hh"

#include <iosfwd>
#include <string>
#include <vector>


//#define DEBUG_SFINDER


/// \brief Build an SV locus graph from read evidence extracted from a contiguous bam segment
///
/// This object is used to build an SV locus graph (SVLocusSet),from a contiguous bam
/// segment, given that sequencing read evidence will be provided with strictly increasing
/// aligned position values from one end to the other of the bam segment in question.
///
/// pos_processor_base:
///
/// This inherits from pos_processor_base to facilitate a "rolling" execution of
/// functions at a defined positional offset less than the position of the most recent
/// read alignment input. These offset functions will (1) trigger the inline graph denoising
/// process (2) clean up buffered read depth data after it is no longer needed.
///
/// Depth Tracking and High-Depth Filtration:
///
/// This object tracks the estimated depth per position summed over all non-tumor samples in
/// the input. This depth estimate is compared to a precomputed expected depth for the chromosome/contig
/// currently being analyzed. At positions where the local read depth is substantially higher than the
/// expected chromosome depth, SV evidence reads are not merged into the SV graph.
///
/// The motivation for skipping high depth regions is that (in an unbiased/non-targeted sequencing assay)
/// high-depth regions indicated reference compressions and other artifacts which typically produce many
/// artifactual SV calls, skipping these regions is an optimization of both FP rate and runtime. Tumor
/// depth is not used because very high depth in the tumor is much more likely to indicate a somatic
/// mutation process, so at least on non-tumor sample is required to infer the location of reference
/// compressions or similarly unreliable regions.
///
struct SVLocusSetFinder : public pos_processor_base
{
    /// This constructs to an immediately usable state following an RAII-like pattern.
    ///
    /// \param scanRegion The genomic region which this SVLocusSetFinder object will translate into an SVLocusGraph
    /// \param bamHeader Bam header containing chromosome details. This is used directly in this object and copied
    ///                  into the SV locus graph.
    SVLocusSetFinder(
        const ESLOptions& opt,
        const GenomeInterval& scanRegion,
        const bam_header_info& bamHeader,
        const reference_contig_segment& refSeq);

    ~SVLocusSetFinder() override
    {
        flush();
    }

    /// \brief Push a new read alignment into the SV graph building process
    ///
    /// \param[in] streamErrorReporter Reference to the error reporter from the stream which produced \p bamRead.
    ///              This is only used to improve the detail of exception messages.
    /// \param[in] bamRead The BAM/CRAM record being pushed into the estimation process.
    /// \param[in] defaultReadGroupIndex The read group index to use in the absence of a BAM read group (RG) tag.
    /// This should effectively be the same as the sample index.
    ///
    void
    update(
        const stream_state_reporter& streamErrorReporter,
        const bam_record& bamRead,
        const unsigned defaultReadGroupIndex);

    /// \brief Provide const access to the SV locus graph that this object is building.
    const SVLocusSet&
    getLocusSet()
    {
        return _svLoci;
    }

    /// \brief Flush any cached values built up during the update process.
    ///
    /// Calling this method should ensure that the SV locus graph reflects all read evidence input so far, and the
    /// graph is represented as compactly as possible.
    void
    flush()
    {
        _stageManager.reset();
    }

    /// \brief Record the time elapsed in the graph building process.
    ///
    /// The SV locus graph object tracks various meta-data related to the graph
    /// including the time spent conducting the graph construction process. This
    /// information is only used for debug/audit and does not impact results.
    ///
    void
    setBuildTime(const CpuTimes& t)
    {
        _svLoci.setBuildTime(t);
    }

private:

    /// Execute logic which is dependent on being a fixed offset from the
    /// the HEAD position in this object's "rolling" positional processing pipeline.
    ///
    /// For this object the HEAD position is defined by the mapping position of the
    /// most recently processed BAM input read.
    ///
    /// For more general background see stage_manager and its associated tests.
    ///
    /// \param stage_no the stage id is used to determine what logic to execute on the given position
    /// \param pos execute stage specific logic on this position number
    void
    process_pos(const int stage_no,
                const pos_t pos) override;

    /// \brief Add the input read to this object's running estimate of read depth per position.
    ///
    /// (see class docs for overview of high depth filtration)
    void
    addToDepthBuffer(
        const bam_record& bamRead);

    enum hack_t
    {
        /// Length in bases on the beginning and the end of scan range which is excluded from in-line graph de-noising
        ///
        /// TODO compute this number from read insert ranges
        REGION_DENOISE_BORDER = 5000
    };

    /////////////////////////////////////////////////
    // data:

    /// Maps the Read Group/Sample Index to a bit indicating whether the index comes from a tumor sample.
    ///
    /// This is required for high-depth filtration because high-depth can only be determined from non-tumor samples.
    const std::vector<bool> _isAlignmentTumor;

    /// The target genome region for this SV locus graph building process
    const GenomeInterval _scanRegion;

    /// A subset of _scanRegion in which the inline graph denoising operation is allowed.
    const GenomeInterval _denoiseRegion;

    /// Helper object used to schedule calls at positions with a defined offset below the current input read's position
    stage_manager _stageManager;

    /// The SV locus graph being built by this object
    SVLocusSet _svLoci;

    /// Track estimated depth per position for the purpose of filtering high-depth regions
    depth_buffer_compressible _positionReadDepthEstimate;

    /// True when the denoising position pointer is within _denoiseRegion
    bool _isInDenoiseRegion;

    /// The graph inline denoising routine has run up to this position, subseqeunt denoising passes should start
    /// from this position up to a chosen window size
    pos_t _denoiseStartPos;

    SVLocusScanner _readScanner;

    /// If true, then track estimated depth/pos and filter out input from very high depth regions
    /// This would typically be true for WGS and false for targeted sequencing.
    bool _isMaxDepthFilter;

    /// Above this depth, read input is not added to the SV Locus graph.
    ///
    /// This depth is supplied from an external estimate
    float _maxDepth;

    const bam_header_info& _bamHeader;
    const reference_contig_segment& _refSeq;
};
