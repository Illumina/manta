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

#pragma once

#include "ESLOptions.hh"

#include "blt_util/bam_record.hh"
#include "blt_util/pos_processor_base.hh"
#include "blt_util/stage_manager.hh"
#include "manta/ReadGroupStatsSet.hh"
#include "manta/SVLocusSet.hh"

#include <iosfwd>
#include <string>
#include <vector>


// estimate an SVLocusSet
//
struct SVLocusSetFinder : public pos_processor_base
{
    SVLocusSetFinder(
            const ESLOptions& opt,
            const GenomeInterval& scanRegion);

    ~SVLocusSetFinder()
    {
        _stageman.reset();
    }

    ///
    /// index is the read group index to use by in the absense of an RG tag
    /// (for now RGs are ignored for the purpose of gathering insert stats)
    ///
    void
    update(const bam_record& read,
           const unsigned defaultReadGroupIndex);

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

private:

    struct CachedReadGroupStats
    {
        CachedReadGroupStats() :
            min(0),
            max(0)
        {}

        double min;
        double max;
    };

    /// this predicate runs any fast tests on the acceptability of a
    /// read for the SVLocus build
    bool
    isReadFiltered(const bam_record& read) const;

    static
    void
    getChimericSVLocus(
        const CachedReadGroupStats& rstats,
        const bam_record& read,
        SVLocus& locus);

    void
    process_pos(const int stage_no,
                const pos_t pos);

    void
    updateDenoiseRegion();

    // TODO -- compute this number from read insert ranges:
    enum hack_t {
        REGION_DENOISE_BORDER = 5000
    };

    /////////////////////////////////////////////////
    // data:
    const ESLOptions& _opt;
    const GenomeInterval _scanRegion;
    GenomeInterval _denoiseRegion;
    stage_manager _stageman;
    ReadGroupStatsSet _rss;
    SVLocusSet _svLoci;

    std::vector<CachedReadGroupStats> _stats;

    bool _isScanStarted;

    bool _isInDenoiseRegion;
    pos_t _denoisePos;
};

