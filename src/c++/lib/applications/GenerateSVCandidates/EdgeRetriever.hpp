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

/// Provide an iterator over edges in a set of SV locus graphs
///
struct EdgeRetriever {
  EdgeRetriever(const SVLocusSet& set, const unsigned graphNodeMaxEdgeCount)
    : _set(set), _graphNodeMaxEdgeCount(graphNodeMaxEdgeCount)
  {
  }

  virtual ~EdgeRetriever() {}

  virtual bool next() = 0;

  const EdgeInfo& getEdge() const { return _edge; }

protected:
  const SVLocusSet& _set;
  const unsigned    _graphNodeMaxEdgeCount;
  EdgeInfo          _edge;
};
