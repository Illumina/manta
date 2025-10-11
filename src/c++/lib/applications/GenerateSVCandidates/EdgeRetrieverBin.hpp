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
///

#pragma once

#include "EdgeRetriever.hpp"

/// Provide an iterator over edges in a set of SV locus graphs
///
/// This facilitates parallelization of graph processing by dividing the graph edges into 'bins',
/// where the edges in each bin have a similar total edge observation counts
///
struct EdgeRetrieverBin final : public EdgeRetriever {
  /// \param[in] graphNodeMaxEdgeCount Filtration parameter for skipping edges from highly connected nodes
  /// (set to zero to disable)
  ///
  /// \param[in] binCount Total number of parallel bins, must be 1 or greater
  ///
  /// \param[in] binIndex Parallel bin id, must be less than binCount
  EdgeRetrieverBin(
      const SVLocusSet& set,
      const unsigned    graphNodeMaxEdgeCount,
      const unsigned    binCount,
      const unsigned    binIndex);

  bool next() override;

private:
  /// Advance the EdgeRetriever::_edge pointer to the first unfiltered edge such that the sum of evidence from
  /// the pointer edge and all previous edges is greater than _beginCount. _headCount will be updated to
  /// reflect the above sum
  ///
  /// In the case that all edges where _headCount is greater than _beginCount are filtered, the function will
  /// set _headCount to _endCount+1 and _edge to an undefined state.
  ///
  void jumpToFirstEdge();

  void advanceEdge();

  /// Provide the observation range for the bin we're retrieving, these values should be constant following
  /// the ctor
  unsigned long _beginCount;
  unsigned long _endCount;

  /// Tracking index of cumulative observation count as we step through the graph
  unsigned long _headCount;
};
