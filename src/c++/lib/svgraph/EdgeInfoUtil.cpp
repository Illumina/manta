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

#include "svgraph/EdgeInfoUtil.hpp"

bool testIsolatedEdge(const SVLocusSet& cset, const EdgeInfo& edge)
{
  if (edge.nodeIndex1 != edge.nodeIndex2) return false;

  const SVLocus& locus(cset.getLocus(edge.locusIndex));
#if 0
    // simple criteria -- make sure there are no other nodes in locus
    return (locus.size() == 1);
#endif

  // search to check to see if any bidirectional edges extend from this node (other than the self-edge):
  const SVLocusNode&       node1(locus.getNode(edge.nodeIndex1));
  const SVLocusEdgeManager node1Manager(node1.getEdgeManager());

  typedef SVLocusEdgesType::const_iterator edgeiter_t;
  edgeiter_t                               edgeIter(node1Manager.getMap().begin());
  const edgeiter_t                         edgeIterEnd(node1Manager.getMap().end());

  EdgeInfo testEdge(edge);

  unsigned edgeCount(0);
  unsigned biEdgeCount(0);
  for (; edgeIter != edgeIterEnd; ++edgeIter) {
    testEdge.nodeIndex2 = (edgeIter->first);
    if (testEdge.nodeIndex1 == testEdge.nodeIndex2) continue;
    edgeCount++;
    if (isBidirectionalEdge(cset, testEdge)) biEdgeCount++;
  }

  const bool isLowBiEdge((biEdgeCount >= 1) && (biEdgeCount <= 2));
  const bool isLowTotalEdge(edgeCount <= 4);
  return (!(isLowBiEdge && isLowTotalEdge));
}
