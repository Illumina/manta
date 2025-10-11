//
// Manta - Structural Variant and Indel Caller
// Copyright (c) 2013-2025 Illumina, Inc.
//
// This program is licensed under the terms of the Polyform strict license
//
// ***As far as the law allows, the software comes as is, without
// any warranty or condition, and the licensor will not be liable
// to you for any damages arising out of these terms or the use
// or nature of the software, under any kind of legal claim.***
//
// You should have received a copy of the PolyForm Strict License 1.0.0
// along with this program.  If not, see <https://polyformproject.org/licenses/strict/1.0.0>.
//
//

/// \file
/// \author Chris Saunders
/// \author Trevor Ramsay
///

#pragma once

#include <memory>

#include "boost/noncopyable.hpp"

#include "blt_util/pos_processor_base.hpp"
#include "blt_util/stage_manager.hpp"
#include "svgraph/GenomeIntervalUtil.hpp"

struct depth_buffer_compressible;
struct SVLocusSet;

/// \brief This objects helps SVLocusSetFinder manage data which needs to be updated as a funciton of the
/// current genome scanning location.
///
/// SVLocusSetFinder is designed to be used by scanning regions of the genome left ot right. Certain data
/// structures in SVLocusSetFinder need to be updated as a function of where that scan process currently is
/// within the target region. These region dependent management tasks are isolated to this ActiveRegionManager
/// object.
///
/// pos_processor_base:
///
/// This inherits from pos_processor_base to facilitate a "rolling" execution of functions at a defined
/// positional offset less than the position of the most recent read alignment input. These offset functions
/// will (1) trigger the inline graph denoising process (2) clean up buffered read depth data after it is no
/// longer needed.
///
struct SVLocusSetFinderActiveRegionManager : public pos_processor_base, private boost::noncopyable {
  /// \param[in] scanRegion The genomic region which the SVLocusSetFinder object will translate into an
  /// SVLocusGraph
  ///
  /// \param[in] svLociPtr Pointer to the SV Locus graph requiring region-based management
  /// \param[in] positionReadDepthEstimatePtr Pointer to depth estimate tracker, which needs to be cleared as
  /// a function of the current active region
  ///
  /// \param[in] denoiseRegionProtectedBorderSize Length in bases of the region protected from denoising. This
  /// region extends back from the highest value provided so far to handle_new_pos_value AND is applied on the
  /// beginning and the end of the scan range, unless the scan range is determined to intersect the edge of a
  /// chromosome
  ///
  SVLocusSetFinderActiveRegionManager(
      const GenomeInterval&                      scanRegion,
      std::shared_ptr<SVLocusSet>                svLociPtr,
      std::shared_ptr<depth_buffer_compressible> positionReadDepthEstimatePtr,
      unsigned                                   denoiseRegionProtectedBorderSize = 5000);

  ~SVLocusSetFinderActiveRegionManager() override { flush(); }

  /// \brief Exposed method from stage_manager. See stage_manager::handle_new_pos_value
  void handle_new_pos_value(const pos_t pos) { _stageManager.handle_new_pos_value(pos); }

  /// \brief Flush any cached values built up during the update process.
  ///
  /// Calling this method should ensure that the SV locus graph reflects all read evidence input so far, and
  /// the graph is represented as compactly as possible.
  void flush() { _stageManager.reset(); }

private:
  /// Execute logic which is dependent on being a fixed offset from the
  /// the HEAD position in this object's "rolling" positional processing pipeline.
  ///
  /// For this object the HEAD position is defined by the mapping position of the
  /// most recently processed BAM input read.
  ///
  /// For more general background see stage_manager and its associated tests.
  ///
  /// \param stage_no the stage id is used to determine what logic to execute on the given position
  ///
  /// \param pos execute stage specific logic on this position number
  void process_pos(const int stage_no, const pos_t pos) override;

  SVLocusSet& _getLocusSet() { return (*_svLociPtr); }

  /////////////////////////////////////////////////
  // data:

  /// The target genome region for this SV locus graph building process
  const GenomeInterval _scanRegion;

  /// Pointer to the SV locus graph object for which region-based management is required
  std::shared_ptr<SVLocusSet> _svLociPtr;

  /// Track estimated depth per position for the purpose of filtering high-depth regions
  std::shared_ptr<depth_buffer_compressible> _positionReadDepthEstimatePtr;

  /// A subset of _scanRegion in which the inline graph denoising operation is allowed.
  const GenomeInterval _denoiseRegion;

  /// Helper object used to schedule calls at positions with a defined offset below the current input read's
  /// position
  stage_manager _stageManager;

  /// True when the denoising position pointer is within _denoiseRegion
  bool _isInDenoiseRegion;

  /// The graph inline denoising routine has run up to this position, subseqeunt denoising passes should start
  /// from this position up to a chosen window size
  pos_t _denoiseStartPos;
};
