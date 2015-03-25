// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Manta
// Copyright (c) 2013-2015 Illumina, Inc.
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

#pragma once

#include "ESLOptions.hh"

#include "blt_util/depth_buffer.hh"
#include "blt_util/pos_processor_base.hh"
#include "blt_util/stage_manager.hh"
#include "htsapi/bam_record.hh"
#include "manta/SVLocusScanner.hh"
#include "svgraph/SVLocusSet.hh"
#include "truth/TruthTracker.hh"

#include <iosfwd>
#include <string>
#include <vector>


//#define DEBUG_SFINDER



/// estimate an SVLocusSet
///
struct SVLocusSetFinder : public pos_processor_base
{
    SVLocusSetFinder(
        const ESLOptions& opt,
        const GenomeInterval& scanRegion,
        const bam_header_info& bamHeader);

    ~SVLocusSetFinder()
    {
        flush();
    }

    /// index is the read group index to use by in the absence of an RG tag
    /// (for now RGs are ignored for the purpose of gathering insert stats)
    ///
    void
    update(
        const bam_record& bamRead,
        const unsigned defaultReadGroupIndex,
        const bam_header_info& bamHeader,
        const reference_contig_segment& refSeq,
        TruthTracker& truthTracker);

    const SVLocusSet&
    getLocusSet()
    {
        return _svLoci;
    }

    void
    setBamHeader(const bam_header_t& header)
    {
        assert(! _isScanStarted);
        _svLoci.header = bam_header_info(header);
        updateDenoiseRegion();
    }

    // flush any cached values built up during the update process
    void
    flush()
    {
        _stageman.reset();
    }

    void
    setBuildTime(const CpuTimes& t)
    {
        _svLoci.setBuildTime(t);
    }

private:

    void
    process_pos(const int stage_no,
                const pos_t pos);

    void
    updateDenoiseRegion();

    void
    addToDepthBuffer(
        const unsigned defaultReadGroupIndex,
        const bam_record& bamRead);

    // TODO -- compute this number from read insert ranges:
    enum hack_t
    {
        REGION_DENOISE_BORDER = 5000    ///< length in bases on the beginning and the end of scan range which is excluded from in-line graph de-noising
    };

    /////////////////////////////////////////////////
    // data:
    const std::vector<bool> _isAlignmentTumor;
    const GenomeInterval _scanRegion;
    GenomeInterval _denoiseRegion;
    stage_manager _stageman;
    SVLocusSet _svLoci;

    depth_buffer_compressible _depth; ///< track depth for the purpose of filtering high-depth regions

    bool _isScanStarted;

    bool _isInDenoiseRegion;
    pos_t _denoisePos;

    SVLocusScanner _readScanner;

    bool _isMaxDepth;
    float _maxDepth;
};

