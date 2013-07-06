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

#include "blt_util/bam_record.hh"
#include "manta/ReadGroupStatsSet.hh"
#include "manta/SVCandidate.hh"
#include "manta/SVLocus.hh"
#include "options/ReadScannerOptions.hh"

#include <string>
#include <vector>


/// consolidate functions which process a read to determine its
/// SV evidence value
///
/// In manta, evidence is scanned (at least) twice: once for SVLocus Graph generation
/// and then once again during hygen/scoring. We need to make sure both of these steps
/// are using the same logic to process read pairs into SV evidence. This class is
/// responsible for the shared logic
///
struct SVLocusScanner
{
    SVLocusScanner(
        const ReadScannerOptions& opt,
        const std::string& statsFilename,
        const std::vector<std::string>& alignmentFilename);

    /// this predicate runs any fast tests on the acceptability of a
    /// read for the SVLocus build
    bool
    isReadFiltered(const bam_record& read) const;

    /// if read supports a chimera candidate return this as a single observation SVLocus object,
    /// else return an empty object.
    ///
    /// \param defaultReadGroupIndex the read group index to use by in the absence of an RG tag
    /// (for now RGs are ignored for the purpose of gathering insert stats)
    ///
    void
    getChimericSVLocus(const bam_record& read,
                       const unsigned defaultReadGroupIndex,
                       SVLocus& locus) const;

    /// get local and remote breakends from read pair
    ///
    /// if remote read is not available, set to NULL and best estimate will be generated
    ///
    void
    getBreakendPair(const bam_record& localRead,
                    const bam_record* remoteReadPtr,
                    const unsigned defaultReadGroupIndex,
                    SVBreakend& localBreakend,
                    SVBreakend& remoteBreakend) const;

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


    static
    void
    getReadBreakendsImpl(
        const CachedReadGroupStats& rstats,
        const bam_record& localRead,
        const bam_record* remoteReadPtr,
        SVBreakend& localBreakend,
        SVBreakend& remoteBreakend,
        known_pos_range2& evidenceRange);

    static
    void
    getChimericSVLocusImpl(
        const CachedReadGroupStats& rstats,
        const bam_record& read,
        SVLocus& locus);

    /////////////////////////////////////////////////
    // data:
    const ReadScannerOptions _opt;
    ReadGroupStatsSet _rss;

    std::vector<CachedReadGroupStats> _stats;
};

