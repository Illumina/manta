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

#include "EdgeOptions.hpp"
#include "EdgeRetriever.hpp"

/// provide an iterator over edges in a set of SV locus graphs
///
/// designed to allow parallelization of the graph processing by
/// dividing iteration into a set of bins with similar total edge
/// observation counts
///
struct EdgeRetrieverLocus final : public EdgeRetriever {
  /// \param locusIndex iterate over all edges of a specific locus
  EdgeRetrieverLocus(
      const SVLocusSet& set, const unsigned graphNodeMaxEdgeCount, const LocusEdgeOptions& opt);

  bool next() override;

private:
  void advanceEdge();

  LocusEdgeOptions _opt;
  bool             _isInit;
};
