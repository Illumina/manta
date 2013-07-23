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
#include "svgraph/SVLocus.hh"
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
    isReadFiltered(const bam_record& bamRead) const;

    /// custom version of proper pair bit test:
    bool
    isProperPair(
        const bam_record& bamRead,
        const unsigned defaultReadGroupIndex) const;

    /// return zero to many SVLocus objects if the read supports any structural variant(s) (detectable by Manta)
    ///
    void
    getSVLoci(
        const bam_record& bamRead,
        const unsigned defaultReadGroupIndex,
        std::vector<SVLocus>& loci) const;

    /// get local and remote breakends for each SV Candidate which can be extracted from a read pair
    ///
    /// if remote read is not available, set remoteReadPtr to NULL and a best estimate will be generated for the remote breakend
    ///
    /// for all candidates, if one breakend is estimated from localRead and one is estimated from remoteRead, then
    /// the local breakend will be placed in candidate bp1 and the remote breakend will be placed in candidate.bp2
    ///
    void
    getBreakendPair(
        const bam_record& localRead,
        const bam_record* remoteReadPtr,
        const unsigned defaultReadGroupIndex,
        std::vector<SVCandidate>& candidates) const;

private:

    struct Range
    {
        Range() :
            min(0),
            max(0)
        {}

        double min;
        double max;
    };

    struct CachedReadGroupStats
    {
        Range breakendRegion;
        Range properPair;
    };


    static
    void
    getReadBreakendsImpl(
        const ReadScannerOptions& opt,
        const CachedReadGroupStats& rstats,
        const bam_record& localRead,
        const bam_record* remoteReadPtr,
        std::vector<SVCandidate>& candidates,
        known_pos_range2& localEvidenceRange);

    static
    void
    getSVLociImpl(
        const ReadScannerOptions& opt,
        const CachedReadGroupStats& rstats,
        const bam_record& bamRead,
        std::vector<SVLocus>& loci);

    /////////////////////////////////////////////////
    // data:
    const ReadScannerOptions _opt;
    ReadGroupStatsSet _rss;

    std::vector<CachedReadGroupStats> _stats;
};

