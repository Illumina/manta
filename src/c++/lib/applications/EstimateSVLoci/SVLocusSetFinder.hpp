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
///

#pragma once

#include "ESLOptions.hpp"
#include "SVLocusSetFinderActiveRegionManager.hpp"

#include "blt_util/depth_buffer.hpp"
#include "htsapi/bam_record.hpp"
#include "htsapi/bam_streamer.hpp"
#include "manta/SVLocusScanner.hpp"
#include "svgraph/SVLocusSet.hpp"

#include "boost/noncopyable.hpp"

#include <iosfwd>
#include <memory>
#include <string>
#include <vector>

//#define DEBUG_SFINDER

/// \brief Build into an SV locus graph from read evidence extracted from a contiguous bam segment
///
/// This object is used to build more information into a given SV locus graph (SVLocusSet),from a contiguous
/// bam segment, given that sequencing read evidence will be provided with strictly increasing aligned
/// position values from one end to the other of the bam segment in question.
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
struct SVLocusSetFinder : private boost::noncopyable {
  /// This constructs to an immediately usable state following an RAII-like pattern.
  ///
  /// \param[in] scanRegion The genomic region which this SVLocusSetFinder object will translate into an
  /// SVLocusGraph
  /// \param[in] refSeqPtr Pointer reference segment corresponding to scanRegion.
  /// \param[in,out] svLociPtr Pointer to SVLocusSet into which all results wil be merged.
  SVLocusSetFinder(
      const ESLOptions&                               opt,
      const GenomeInterval&                           scanRegion,
      const std::shared_ptr<reference_contig_segment> refSeqPtr,
      std::shared_ptr<SVLocusSet>                     svLociPtr);

  /// \brief Push a new read alignment into the SV graph building process
  ///
  /// \param[in] streamErrorReporter Reference to the error reporter from the stream which produced
  /// \p bamRead. This is only used to improve the detail of exception messages.
  /// \param[in] bamRead The BAM/CRAM record being pushed into the estimation process.
  /// \param[in] defaultReadGroupIndex The read group index to use in the absence of a BAM read group (RG)
  /// tag. This should effectively be the same as the sample index.
  ///
  void update(
      const stream_state_reporter& streamErrorReporter,
      const bam_record&            bamRead,
      const unsigned               defaultReadGroupIndex);

  /// \brief Provide const access to the SV locus graph that this object is building.
  const SVLocusSet& getLocusSet() const { return *_svLociPtr; }

  /// \brief Flush any cached values built up during the update process.
  ///
  /// Calling this method should ensure that the SV locus graph reflects all read evidence input so far, and
  /// the graph is represented as compactly as possible.
  void flush() { _regionManager.flush(); }

private:
  /// \brief Add the input read to this object's running estimate of read depth per position.
  ///
  /// (see class docs for overview of high depth filtration)
  void addToDepthBuffer(const bam_record& bamRead);

  const bam_header_info& _bamHeader() const { return getLocusSet().getBamHeader(); }

  const reference_contig_segment& _refSeq() const { return *(_refSeqPtr.get()); }

  SVLocusSet& _getLocusSet() { return *_svLociPtr; }

  /////////////////////////////////////////////////
  // data:

  /// Maps the Read Group/Sample Index to a bit indicating whether the index comes from a tumor sample.
  ///
  /// This is required for high-depth filtration because high-depth can only be determined from non-tumor
  /// samples.
  const std::vector<bool> _isAlignmentTumor;

  /// The target genome region for this SV locus graph building process
  const GenomeInterval _scanRegion;

  const std::shared_ptr<reference_contig_segment> _refSeqPtr;

  /// Pointer to the SV locus graph being built into by this object
  std::shared_ptr<SVLocusSet> _svLociPtr;

  /// Track estimated depth per position for the purpose of filtering high-depth regions
  std::shared_ptr<depth_buffer_compressible> _positionReadDepthEstimatePtr;

  /// Region manager handles all region-dependent data triggers (denoising, cleaning unused depth positions,
  /// etc..)
  SVLocusSetFinderActiveRegionManager _regionManager;

  SVLocusScanner _readScanner;

  /// If true, then track estimated depth/pos and filter out input from very high depth regions
  /// This would typically be true for WGS and false for targeted sequencing.
  bool _isMaxDepthFilter;

  /// Above this depth, read input is not added to the SV Locus graph.
  ///
  /// This depth is supplied from an external estimate
  float _maxDepth;
};
