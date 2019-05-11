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
