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
// <https://github.com/sequencing/licenses/>
//

///
/// \author Chris Saunders
///

#pragma once

#include "ESLOptions.hh"

#include "blt_util/bam_record.hh"
#include "blt_util/pos_processor_base.hh"
#include "blt_util/stage_manager.hh"
#include "manta/SVLocusScanner.hh"
#include "svgraph/SVLocusSet.hh"

#include <iosfwd>
#include <string>
#include <vector>


//#define DEBUG_SFINDER



// estimate an SVLocusSet
//
struct SVLocusSetFinder : public pos_processor_base
{
    SVLocusSetFinder(
        const ESLOptions& opt,
        const GenomeInterval& scanRegion);

    ~SVLocusSetFinder()
    {
        flush();
    }

    ///
    /// index is the read group index to use by in the absence of an RG tag
    /// (for now RGs are ignored for the purpose of gathering insert stats)
    ///
    void
    update(const bam_record& bamRead,
           const unsigned defaultReadGroupIndex,
           const std::string& bkptRef,
           const std::map<std::string, int32_t>& chromToIndex);

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
        _svLoci.addAnomCount(_anomCount);
        _svLoci.addNonAnomCount(_nonAnomCount);
        _stageman.reset();

        _anomCount=0;
        _nonAnomCount=0;
    }


private:

    void
    process_pos(const int stage_no,
                const pos_t pos);

    void
    updateDenoiseRegion();

    // TODO -- compute this number from read insert ranges:
    enum hack_t
    {
        REGION_DENOISE_BORDER = 5000    ///< length in bases on the begining and the end of scan range which is excluded from in-line graph de-noising
    };

    /////////////////////////////////////////////////
    // data:
    const GenomeInterval _scanRegion;
    GenomeInterval _denoiseRegion;
    stage_manager _stageman;
    SVLocusSet _svLoci;

    bool _isScanStarted;

    bool _isInDenoiseRegion;
    pos_t _denoisePos;

    SVLocusScanner _readScanner;

    unsigned _anomCount;
    unsigned _nonAnomCount;
};

