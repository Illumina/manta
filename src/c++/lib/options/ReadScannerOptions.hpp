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

#pragma once

struct ReadScannerOptions {
  ReadScannerOptions() {}

  /// \brief Reads with MAPQ below this value are filtered out during SV locus generation and many subsequent
  /// candidate creation and scoring steps
  ///
  unsigned minMapq = 15;

  /// \brief Reads with MAPQ below this value may be considered in the normal sample as part of the 'tier2'
  /// evidence to help disprove a somatic candidate via the existence of weak candidate evidence in the normal
  unsigned minTier2Mapq = 5;

  /// \brief Probability used to create fragment size quantile range for reported breakend regions
  float breakendEdgeQuantileProb = 0.25f;

  /// \brief Probability used to create fragment size quantile range for reported breakend regions, this
  /// version is reserved for large-scale SV candidates
  float largeScaleEventBreakendEdgeQuantileProb = 0.1f;

  /// \brief Probability used to create fragment size quantile range for anomalous read pair detection during
  /// SV discovery
  ///
  /// Treat a paired read as non-anomalous (ie. "proper") during SV discovery if orientation is non-anomalous
  /// and implied fragment size is within the [x, 1-x], x=properPairQuantileProb quantile range of the
  /// fragment length distribution.
  float properPairQuantileProb = 0.01f;

  /// \brief Probability used to create fragment size quantile range for supporting read pair evidence during
  /// SV scoring
  ///
  /// Add a read pair to an SV's pool of evidence to evaluate for paired-read support during the scoring phase
  /// if the orientation is non-anomalous and the implied fragment size is within the [x, 1-x],
  /// x=evidenceQuantileProb quantile range of the fragment length
  float evidenceTrimQuantileProb = 0.15f;

  /// Shadow read support is searched upstream of a breakend to a distance of X*Y, where
  /// X=quantile(1-shadowSearchDistanceQuantileProb) on the fragment length distribution
  /// Y=shadowSearchDistanceFactor
  float shadowSearchDistanceQuantileProb = 0.05f;

  /// Shadow read support is searched upstream of a breakend to a distance of X*Y, where
  /// X=quantile(shadowSearchDistanceProb) on the fragment length distribution
  /// Y=shadowSearchDistanceFactor
  float shadowSearchDistanceFactor = 1.2f;

  /// \brief Ignore indels smaller than this when building graph, constructing candidates and scoring output
  unsigned minCandidateVariantSize = 10;

  /// if minCandidateVariantSize is set higher than this value, then we ignore non-specific assembly evidence
  /// (like soft-clip and poor alignment) during candidate generation
  unsigned maxCandidateSizeForLocalAssmEvidence = 100;

  /// whenever a breakend is predicted from a read pair junction, the predicted breakend
  /// range should be no smaller than this:
  unsigned minPairBreakendSize = 40;

  /// \brief Minimum length required of a poorly aligned read end-segment for it to be treated as SV evidence.
  unsigned minSemiAlignedMismatchLen = 8;

  /// \brief Minimum length of a cis SV candidate in RNA data
  unsigned minRNACisLength = 100000;
  /// \brief Minimum length of any SV candidate in RNA data
  unsigned minRNALength = 1000;

  /// \brief Minimum mapping quality for shadow mate used to build SV adjacency graph
  unsigned minSingletonMapqGraph = 30;

  /// \brief Minimum mapping quality for shadow mate used for candidate assembly and scoring
  unsigned minSingletonMapqCandidates = 15;

  // \brief If true, consider an overlapping read pair as evidence.
  bool useOverlapPairEvidence = false;

  /// \brief If true, do not treat reads with the 'proper pair' bit set as SV evidence.
  ///
  /// This is typically set true for RNA-Seq analysis, where proper-pair is used to signal intron-spanning
  /// pairs.
  bool isIgnoreAnomProperPair = false;

  /// \brief The maximum depth factor at which input reads are considered in graph creation/assembly, etc.
  ///
  /// 'depth factor' is the locus' multiple of the expected chromosome depth.
  float maxDepthFactor = 12;

  /// \brief The maximum depth factor over a whole locus. If any site in the locus exceeds this depth, remote
  /// read
  /// retrieval is not used (ie. MAPQ0 chimera mates retrieved for large insertion assembly)
  ///
  /// 'depth factor' is the locus' multiple of the expected chromosome depth.
  float maxLocalDepthFactorForRemoteReadRetrieval = 7;
};
