//
// Manta - Structural Variant and Indel Caller
// Copyright (c) 2013-2019 Illumina, Inc.
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
/// \author Ole Schulz-Trieglaff
/// \author Bret Barnes
///

#pragma once

#include <string>
#include <vector>

#include "blt_util/LinearScaler.hpp"
#include "htsapi/align_path_bam_util.hpp"
#include "htsapi/bam_header_info.hpp"
#include "htsapi/bam_record.hpp"
#include "htsapi/bam_record_util.hpp"
#include "htsapi/bam_streamer.hpp"
#include "manta/ReadGroupStatsSet.hpp"
#include "manta/SVCandidate.hpp"
#include "manta/SVLocusEvidenceCount.hpp"
#include "options/ReadScannerOptions.hpp"
#include "svgraph/SVLocus.hpp"
#include "svgraph/SVLocusSampleCounts.hpp"

namespace FragmentSizeType {
/// \brief Discrete categories used to label DNA fragment length relative to the expected fragment size
/// distribution.
enum index_t {
  COMPRESSED,  ///< The fragment is anomalously small
  NORMAL,      ///< The fragment is non-anomalous
  UNKNOWN,     ///< The fragment size category is unknown
  CLOSE,   ///< The fragment is anomalous, but close enough to the non-anomalous threshold that it is treated
           ///< as possibly non-anomalous.
  DISTANT  ///< The fragment is clearly anomalous
};
}  // namespace FragmentSizeType

/// \brief This object centralizes the weights applied to different sources of SV evidence to indicate
/// relative confidence in the accuracy of each observation.
///
/// The evidence weights used in the SVLocus Graph represent an abstract strength of evidence supporting each
/// edge. The evidence weights are reduced to small counts so that the weights can be used as part of a
/// low-memory graph representation. The weights can be arbitrarily scaled collectively without changing the
/// result, only the  relative values are important.
///
/// Evidence weights express the confidence in an evidence source for discovery and graph building, when SV's
/// are actually scored a much more detailed process is used to derive confidence in the SV based on observed
/// evidence.
struct SVObservationWeights {
  /// Evidence weight for one average-case observation. This acts as a baseline scale for the whole weighting
  /// scheme.
  static const unsigned observation = 3;

  // input evidence:
  static const unsigned readPair          = observation;  ///< Weight used for anomalous read pairs
  static const unsigned closeReadPair     = 1;
  static const unsigned internalReadEvent = observation;  ///< Weight used for indels, soft-clip, etc.
};

struct ReadScannerDerivOptions {
  ReadScannerDerivOptions(const ReadScannerOptions& opt, const bool initIsTranscriptStrandKnown)
    : isSmallCandidates(opt.minCandidateVariantSize <= opt.maxCandidateSizeForLocalAssmEvidence),
      beforeBreakend(opt.minPairBreakendSize / 2),
      afterBreakend(opt.minPairBreakendSize - beforeBreakend),
      isTranscriptStrandKnown(initIsTranscriptStrandKnown)
  {
  }

  const bool  isSmallCandidates;
  const pos_t beforeBreakend;
  const pos_t afterBreakend;

  /// \brief True if running in RNA-Seq mode with a stranded input RNA assay
  const bool isTranscriptStrandKnown;
};

/// \brief Consolidate functions which test aligned reads for SV evidence
///
/// In manta, evidence is scanned (at least) twice: once for SVLocus Graph generation
/// and then once again during hygen/scoring. We need to make sure both of these steps
/// are using the same logic to process read pairs into SV evidence. This class is
/// responsible for the shared logic
///
struct SVLocusScanner {
  SVLocusScanner(
      const ReadScannerOptions&       opt,
      const std::string&              statsFilename,
      const std::vector<std::string>& alignmentFilename,
      const bool                      isTranscriptStrandKnown = false);

  /// QC check if the read length implied by cigar string matches the length of read sequence
  ///
  /// \param[in] alignmentStream This is used (only) to improve the detail of error messages.
  /// \param bamRead
  static void checkReadSize(const stream_state_reporter& alignmentStream, const bam_record& bamRead);

  unsigned getMinMapQ() const { return _opt.minMapq; }

  unsigned getMinTier2MapQ() const { return _opt.minTier2Mapq; }

  /// \brief Test \p bamRead to determine if it is from a paired-read DNA fragment with an anomalous
  /// orientation or inferred fragment size
  ///
  /// Note that in the general case, Manta disregards the bam spec's 'proper pair' flag on input alignments
  /// and instead directly determines anomalous state.
  bool isAnomalousReadPair(const bam_record& bamRead, const unsigned defaultReadGroupIndex) const;

  /// \breif Test \p bamRead to determine if it is from a paired-end DNA fragment with anomalous orientation
  /// or size, excluding the case of an anomalously short standard orientation read pair
  ///
  /// According to this method only mapped read pairs can be anomalous. All single read anomalies (SA tags,
  /// CIGAR, semi-aligned) are detected elsewhere.
  bool isNonCompressedAnomalousReadPair(
      const bam_record& bamRead, const unsigned defaultReadGroupIndex) const;

  /// \brief Test \p bamAlign for large indel segments
  ///
  /// \return True if \p bamAlign contains any indel segments above a critical size threshold
  bool isLocalIndelEvidence(const SimpleAlignment& bamAlign) const;

  /// true for reads with semi-aligned and soft-clipped edges
  ///
  /// "semi-aligned" here means the edge has a high density of mismatches
  bool isSemiAlignedEvidence(
      const bam_record&               bamRead,
      const SimpleAlignment&          bamAlign,
      const reference_contig_segment& refSeq) const;

  /// \brief Return true if the read indicates the presence of a small SV?
  ///
  /// this function flags reads which could contribute to a local small-variant assembly
  /// but would not otherwise be caught by the proper pair function
  ///
  /// "small" here is relative -- it means any event at a size where read pair evidence will not be dominant
  ///
  /// Note that the thresholds in this function are more stringent than the equivalent scan used to
  /// pick up reads prior to assembly -- in this case false positives could clog up the graph and
  /// interfere with larger event discovery if not kept under control
  bool isLocalAssemblyEvidence(const bam_record& bamRead, const reference_contig_segment& refSeq) const;

  /// \brief A fast test to eliminate reads which are very unlikely to contribute any SV or indel evidence
  ///
  /// \param[in,out] incountsPtr If not nullptr, the pointed object is appended with statistics on the SV
  /// evidence type associated with \p bamRead
  ///
  /// \return True if it is probable that \p bamRead provides evidence for an SV or indel
  bool isSVEvidence(
      const bam_record&               bamRead,
      const unsigned                  defaultReadGroupIndex,
      const reference_contig_segment& refSeq,
      SVLocusEvidenceCount*           incountsPtr = nullptr) const;

  /// return zero to many SVLocus objects if the read supports any
  /// structural variant(s) (detectable by manta)
  ///
  /// \param defaultReadGroupIndex the read group index to use in the absence of an RG tag
  /// (for now RGs are ignored for the purpose of gathering insert stats)
  ///
  void getSVLoci(
      const bam_record&               bamRead,
      const unsigned                  defaultReadGroupIndex,
      const bam_header_info&          bamHeader,
      const reference_contig_segment& refSeq,
      std::vector<SVLocus>&           loci,
      SampleEvidenceCounts&           eCounts) const;

  /// \brief Find all candidate SV observations from a single input read pair
  ///
  /// If the alignment record of the remote read from the read pair is not available, set remoteReadPtr to
  /// nullptr and a best estimate will be generated for the remote breakend.
  ///
  /// For all SVObservations generated, if one breakend is estimated from localRead and one is estimated from
  /// remoteRead, then the local breakend will be placed in candidate bp1 and the remote breakend will be
  /// placed in candidate.bp2
  ///
  /// \param[out] candidates Report all candidates found in the input read pair to this structure. This vector
  /// is cleared on input.
  void getBreakendPair(
      const bam_record&               localRead,
      const bam_record*               remoteReadPtr,
      const unsigned                  defaultReadGroupIndex,
      const bam_header_info&          bamHeader,
      const reference_contig_segment& localRefSeq,
      const reference_contig_segment* remoteRefSeqPtr,
      std::vector<SVObservation>&     candidates) const;

  /// \brief Get the distance upstream of a breakend in which shadow read support will be searched for
  ///
  /// \TODO if read groups are decoupled from samples, then this value needs to be revisited to reflect all
  /// read groups in the sample.
  int getShadowSearchDistance(const unsigned defaultReadGroupIndex) const
  {
    return _stats[defaultReadGroupIndex].shadowSearchDistance;
  }

  /// \brief Provide read-only access to the fragment length distribution for each read-group/sample
  ///
  const SizeDistribution& getFragSizeDistro(const unsigned defaultReadGroupIndex) const
  {
    return _rss.getStats(defaultReadGroupIndex).fragStats;
  }

  struct Range {
    double min = 0;
    double max = 0;
  };

  const Range& getEvidencePairRange(const unsigned readGroupIndex) const
  {
    return _stats[readGroupIndex].evidencePair;
  }

  const Range& getExtremeFifthRange() const { return _fifthPerc; }

  /// \brief Summary statistics derived from the full fragment length distribution.
  ///
  /// Instead of transmitting the full fragment distribution for each read group, this object holds the
  /// critical statistics computed from that distribution which are used during discovery and variant calling.
  ///
  /// Note that conceptually, there is one fragment length distribution per "Read Group", but in the current
  /// methods the concept of read group is forced to have a 1-1 relationship with each sample.
  struct CachedReadGroupStats {
    /// \brief Fragment size range used for the purpose of creating SVLocusGraph regions.
    Range breakendRegion;

    /// \brief Fragment size range used for the purpose of creating SVLocusGraph regions,this range is used
    /// exclusively for large scale events (non-deletion or deletion above a threshold size).
    Range largeScaleEventBreakendRegion;

    /// \brief Fragment size range used to determine if a fragment is considered non-anomalous (ie. "proper")
    /// during SV discovery.
    Range properPair;

    /// \brief Fragment size range used to determine if a fragment will be evaluated as paired-read support
    /// during SV scoring.
    Range evidencePair;

    /// \brief Fragment size range over the [0.05,0.95] quantiles of the fragment size distribution.
    Range fifthPerc;

    /// \brief Shadow read support for a breakend is searched for from the breakend location to this distance
    /// upstream.
    int shadowSearchDistance = 0;

    /// \brief Beyond the properPair anomalous threshold, this threshold is used to distinguish 'close' and
    /// 'far' pairs for the purpose of setting the anomalous pair evidence weight
    int minDistantFragmentSize = 0;

    /// \brief Used to set expanded breakend sizes for large events \TODO more detail?
    LinearScaler<int> largeEventRegionScaler;
  };

private:
  /// \brief Classify read pair fragment sizes into a pre-defined size category.
  ///
  /// Note this assumes a mapped read pair. Add a check on this pre-condition if making this method
  /// non-private.
  FragmentSizeType::index_t _getFragmentSizeType(
      const bam_record& bamRead, const unsigned defaultReadGroupIndex) const;

  /// test whether a fragment is significantly larger than expected
  ///
  /// this function is useful to eliminate reads which fail the ProperPair test
  /// but are still very small
  ///
  /// assumes a mapped read pair -- check this in code if making this method non-private
  bool _isLargeFragment(const bam_record& bamRead, const unsigned defaultReadGroupIndex) const;

  /////////////////////////////////////////////////
  // data:
  const ReadScannerOptions      _opt;
  const ReadScannerDerivOptions _dopt;
  ReadGroupStatsSet             _rss;

  std::vector<CachedReadGroupStats> _stats;

  /// extreme 5th-95th percentiles over all read groups:
  Range _fifthPerc;
};
