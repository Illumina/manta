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
/// \author Trevor Ramsay
///

#include "testSVLocusUtil.hpp"

void locusAddPair(
    SVLocus&       locus,
    const int32_t  tid1,
    const int32_t  beginPos1,
    const int32_t  endPos1,
    const int32_t  tid2,
    const int32_t  beginPos2,
    const int32_t  endPos2,
    const bool     bothLocal,
    const unsigned count)
{
  const NodeIndexType nodePtr1 = locus.addNode(GenomeInterval(tid1, beginPos1, endPos1));
  const NodeIndexType nodePtr2 = locus.addNode(GenomeInterval(tid2, beginPos2, endPos2));
  const unsigned      remoteCount(bothLocal ? count : 0);
  locus.linkNodes(nodePtr1, nodePtr2, count, remoteCount);
}
