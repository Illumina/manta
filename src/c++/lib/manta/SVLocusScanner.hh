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
/// \author Ole Schulz-Trieglaff
/// \author Bret Barnes
///

#pragma once

#include "blt_util/bam_record.hh"
#include "blt_util/bam_record_util.hh"
#include "blt_util/LinearScaler.hh"
#include "manta/ReadGroupStatsSet.hh"
#include "manta/SVCandidate.hh"
#include "svgraph/SVLocus.hh"
#include "options/ReadScannerOptions.hh"
#include "truth/TruthTracker.hh"

#include <string>
#include <vector>


namespace FragmentSizeType
{
enum index_t
{
    COMPRESSED,
    NORMAL,
    VERYCLOSE,
    CLOSE,
    DISTANT
};
}


/// The counts in the SVLocus Graph represent an abstract weight of evidence supporting each edge/node.
///
/// To support large and small-scale evidence in a single graph, we need to allow for different weightings
/// for different evidence types
///
struct SVObservationWeights
{
    // noise reduction:
    static const unsigned observation = 3; ///< 'average' observation weight, this is used to scale noise filtration, but not for any evidence type

    // input evidence:
    static const unsigned readPair = observation;
    static const unsigned closeReadPair = 1;
    static const unsigned veryCloseReadPair = 1;
    static const unsigned internalReadEvent = observation; ///< indels, soft-clip, etc.
};



bool
isGoodShadow(
    const bam_record& bamRead,
    const uint8_t lastMapq,
    const std::string& lastQname,
    const double minSingletonMapq);


struct ReadScannerDerivOptions
{
    ReadScannerDerivOptions(const ReadScannerOptions& opt) :
        beforeBreakend(opt.minPairBreakendSize/2),
        afterBreakend(opt.minPairBreakendSize-beforeBreakend)
    {}

    const pos_t beforeBreakend;
    const pos_t afterBreakend;
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

    /// this predicate runs isReadFiltered without the mapq components
    static
    bool
    isReadFilteredCore(const bam_record& bamRead)
    {
        if      (bamRead.is_filter()) return true;
        else if (bamRead.is_dup()) return true;
        else if (bamRead.is_secondary()) return true;
        else if (bamRead.is_supplement()) return true;
        return false;
    }

    /// this predicate runs any fast tests on the acceptability of a
    /// read for the SVLocus build
    /// Tests also for low mapq
    bool
    isReadFiltered(const bam_record& bamRead) const
    {
        if      (isReadFilteredCore(bamRead)) return true;
        else if (bamRead.map_qual() < _opt.minMapq) return true;
        return false;
    }

    unsigned
    getMinMapQ() const
    {
        return _opt.minMapq;
    }

    /// custom version of proper pair bit test:
    bool
    isProperPair(
        const bam_record& bamRead,
        const unsigned defaultReadGroupIndex) const;

    /// fragments sizes get thrown is serveral pre-defined categories:
    FragmentSizeType::index_t
    getFragmentSizeType(
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

    /// return true if the read is anomalous, for any anomaly type besides being a short innie read:
    bool
    isNonCompressedAnomalous(
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
        const bam_record& bamRead,
        const reference_contig_segment& refSeq) const;

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
        const std::map<std::string, int32_t>& chromToIndex,
        const reference_contig_segment& refSeq,
        std::vector<SVLocus>& loci,
        TruthTracker& truthTracker) const;

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
        const std::map<std::string, int32_t>& chromToIndex,
        const reference_contig_segment& localRefSeq,
        const reference_contig_segment* remoteRefSeqPtr,
        std::vector<SVObservation>& candidates,
        TruthTracker& truthTracker) const;

    /// provide direct access to the frag distro for
    /// functions which can't be cached
    ///
    const SizeDistribution&
    getFragSizeDistro(
        const unsigned defaultReadGroupIndex) const
    {
        return _rss.getStats(defaultReadGroupIndex).fragStats;
    }

    struct Range
    {
        Range() :
            min(0),
            max(0)
        {}

        double min;
        double max;
    };

    const Range&
    getEvidencePairRange(const unsigned readGroupIndex) const
    {
        return _stats[readGroupIndex].evidencePair;
    }

    struct CachedReadGroupStats
    {
        CachedReadGroupStats() :
            minDistantFragmentSize(0),
            minCloseFragmentSize(0),
            minVeryCloseFragmentSize(0)
        {}

        /// fragment size range assumed for the purpose of creating SVLocusGraph regions
        Range breakendRegion;

        /// fragment size range assumed for the purpose of creating SVLocusGraph regions,
        /// this range is used exclusively for large scale events (non-deletion or deletion above a threshold size):
        Range largeScaleEventBreakendRegion;

        /// fragment size range used to determine if a read is anomalous
        Range properPair;

        Range evidencePair;

        int minDistantFragmentSize; ///< beyond the properPair anomalous threshold, there is a threshold to distinguish close and far pairs for the purpose of evidence weight
        int minCloseFragmentSize; ///< beyond the properPair anomalous threshold, there is a threshold to distinguish 'really-close' and 'close' pairs for the purpose of evidence weight
        int minVeryCloseFragmentSize;

        //LinearScaler<int> veryCloseEventScaler; ///< used to scale down breakend size as fragments get smaller

        LinearScaler<int> largeEventRegionScaler; ///< used to set expanded breakend sizes for large events
    };

private:

    /////////////////////////////////////////////////
    // data:
    const ReadScannerOptions _opt;
    const ReadScannerDerivOptions _dopt;
    ReadGroupStatsSet _rss;

    std::vector<CachedReadGroupStats> _stats;

//    std::string lastQname;
//    uint8_t lastMapq;
};

