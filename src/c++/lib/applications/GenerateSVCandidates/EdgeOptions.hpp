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

/// options for SVLocusGraph edge iteration and noise edge filtration
struct LocusEdgeOptions {
  /// If isLocusIndex, report this locus only
  unsigned locusIndex = 0;

  /// If true, generate candidates for all edges touching a specific node in one locus. Assumes isLocusIndex
  /// is true.
  bool     isNodeIndex1 = false;
  unsigned nodeIndex1   = 0;

  /// If true, generate candidates for only the edge from node1 to node2 in one locus. Assumes isLocusIndex
  /// and isNodeIndex1 are true.
  bool     isNodeIndex2 = false;
  unsigned nodeIndex2   = 0;
};

/// Options for SVLocusGraph edge iteration and noise edge filtration
struct EdgeOptions {
  unsigned binCount = 1;  ///< divide all edges in the graph into binCount bins of approx equal complexity
  unsigned binIndex = 0;  ///< out of binCount bins, iterate through the edges in this bin only

  /// If true, generate candidates for a specific SVgraph locus only, and ignore binCount/binIndex
  bool             isLocusIndex = false;
  LocusEdgeOptions locusOpt;

  /// If both nodes of an edge have an edge count higher than this, then skip evaluation of this edge, set to
  /// 0 to turn this filtration off
  unsigned graphNodeMaxEdgeCount = 10;
};
