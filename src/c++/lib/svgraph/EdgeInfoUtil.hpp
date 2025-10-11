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

#include "svgraph/EdgeInfo.hpp"
#include "svgraph/SVLocusSet.hpp"

inline bool isBidirectionalEdge(const SVLocusSet& cset, const EdgeInfo& edge)
{
  const unsigned minEdgeCount(cset.getMinMergeEdgeCount());

  const SVLocus& locus(cset.getLocus(edge.locusIndex));

  return (
      (locus.getEdge(edge.nodeIndex1, edge.nodeIndex2).getCount() >= minEdgeCount) &&
      (locus.getEdge(edge.nodeIndex2, edge.nodeIndex1).getCount() >= minEdgeCount));
}

/// determine if this is a self-edge of a node with no other (bidirectional-pass) edges:
bool testIsolatedEdge(const SVLocusSet& cset, const EdgeInfo& edge);
