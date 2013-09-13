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


/// The counts in the SVLocus Graph represent an abstract weight of evidence supporting each edge/node.
///
/// To support large and small-scale evidence in a single graph, we need to allow for different weightings
/// for different evidence types
///
struct SVObservationWeights
{
    // input evidence:
    static const unsigned readPair = 2;
    static const unsigned closeReadPair = 1;
    static const unsigned internalReadEvent = 2; ///< indels, soft-clip, etc.

    static const float closePairFactor; ///< fragments within this factor of the minimum size cutoff are treated as 'close' pairs and receive a modified evidence count

    // noise reduction:
    static const unsigned observation = 2; ///< 'average' observation weight, this is used to scale noise filtration, but not for any evidence type
};




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
    /// Tests also for low mapq
    bool
    isReadFiltered(const bam_record& bamRead) const;

    /// custom version of proper pair bit test:
    bool
    isProperPair(
        const bam_record& bamRead,
        const unsigned defaultReadGroupIndex) const;

    /// test whether a fragment is significantly larger than expected
    ///
    /// this function is useful to eliminate reads which fail the ProperPair test
    /// but are still very small
    ///
    bool
    isLargeFragment(
        const bam_record& bamRead,
        const unsigned defaultReadGroupIndex) const;

    /// \brief is the read likely to indicate the presence of a small SV?
    ///
    /// this function flags reads which could contribute to a local small-variant assembly
    /// but would not otherwise be caught by the proper pair function
    ///
    /// "small" here is relative -- it means any event at a size where read pair evidence will not be dominant
    ///
    /// Note that the thresholds in this function are more stringent than the equivalent scan used to
    /// pick up reads prior to assembly -- in this case false positives could clog up the graph and
    /// interfere with larger event discovery if not kept under control
    bool
    isLocalAssemblyEvidence(
        const bam_record& bamRead) const;


    /// test for semi-alignedness
    /// Dummy implementation, returns true for now
    bool
    isSemiAligned(
        const bam_record& bamRead) const;

    bool
    isClipped(
        const bam_record& bamRead) const;

    unsigned
    getClipLength(
        const bam_record& bamRead) const;


    /// return zero to many SVLocus objects if the read supports any
    /// structural variant(s) (detectable by manta)
    ///
    /// \param defaultReadGroupIndex the read group index to use in the absence of an RG tag
    /// (for now RGs are ignored for the purpose of gathering insert stats)
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
        CachedReadGroupStats() :
            minFarFragmentSize(0)
        {}

        /// fragment size range assumed for the purpose of creating SVLocusGraph regions
        Range breakendRegion;

        /// fragment size range used to determine if a read is anomalous
        Range properPair;

        Range evidencePair;

        int minFarFragmentSize; ///< beyond the properPair anomalous threshold, there is a threshold to distinguish near and far pairs for the purpose of evidence weight
    };

private:

    /////////////////////////////////////////////////
    // data:
    const ReadScannerOptions _opt;
    ReadGroupStatsSet _rss;

    std::vector<CachedReadGroupStats> _stats;
};

