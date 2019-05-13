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

#include "EdgeRetrieverBin.hpp"

#include <cassert>

//#define DEBUG_EDGER

#ifdef DEBUG_EDGER
#include <iostream>
#include "blt_util/log.hpp"
#endif

/// \brief When \p totalCount is subdivided into \p binCount approximately even bins, return the 0-indexed
/// count which starts the zero-indexed \p binIndex bin.
///
static unsigned long getBoundaryCount(const double binCount, const double binIndex, const double totalCount)
{
  return static_cast<unsigned>(std::floor((totalCount * binIndex) / binCount));
}

EdgeRetrieverBin::EdgeRetrieverBin(
    const SVLocusSet& set,
    const unsigned    graphNodeMaxEdgeCount,
    const unsigned    binCount,
    const unsigned    binIndex)
  : EdgeRetriever(set, graphNodeMaxEdgeCount), _headCount(0)
{
  assert(binCount > 0);
  assert(binIndex < binCount);

  const unsigned long totalObservationCount(_set.totalObservationCount());
  _beginCount = (getBoundaryCount(binCount, binIndex, totalObservationCount));
  _endCount   = (getBoundaryCount(binCount, binIndex + 1, totalObservationCount));

#ifdef DEBUG_EDGER
  log_os << "EDGER: binIndex,binCount,beginCount,endCount: " << binIndex << " " << binCount << " "
         << _beginCount << " " << _endCount << "\n";
#endif
}

void EdgeRetrieverBin::jumpToFirstEdge()
{
  typedef SVLocusEdgesType::const_iterator edgeiter_t;

  const bool     isFilterNodes(_graphNodeMaxEdgeCount > 0);
  const unsigned setSize(_set.size());

  bool isLastFiltered(false);

  // first catch headCount up to the begin edge if required:
  while (true) {
    if (isLastFiltered && (_edge.locusIndex == setSize)) {
      _headCount       = (_endCount + 1);
      _edge.locusIndex = 0;
      _edge.nodeIndex1 = 0;
      _edge.nodeIndex2 = 0;
      return;
    }

    assert(_edge.locusIndex < setSize);

    const SVLocus& locus(_set.getLocus(_edge.locusIndex));
    const unsigned locusObservationCount(locus.totalObservationCount());

    if ((_headCount + locusObservationCount) <= _beginCount) {
      // skip over the locus in this case:
      _headCount += locusObservationCount;
    } else {
      // In this case we expect to (usually) find the start edge within the locus
      const unsigned locusSize(locus.size());

      while (_edge.nodeIndex1 < locusSize) {
        const SVLocusNode& node1(locus.getNode(_edge.nodeIndex1));
        const bool         isEdgeFilterNode1(isFilterNodes && (node1.size() > _graphNodeMaxEdgeCount));

        const SVLocusEdgeManager node1Manager(node1.getEdgeManager());
        edgeiter_t               edgeIter(node1Manager.getMap().lower_bound(_edge.nodeIndex1));
        const edgeiter_t         edgeiterEnd(node1Manager.getMap().cend());

        for (; edgeIter != edgeiterEnd; ++edgeIter) {
          unsigned   edgeCount(edgeIter->second.getCount());
          const bool isSelfEdge(edgeIter->first == _edge.nodeIndex1);
          if (!isSelfEdge) edgeCount += locus.getEdge(edgeIter->first, _edge.nodeIndex1).getCount();

          _headCount += edgeCount;
          isLastFiltered = false;

          if (_headCount > _beginCount) {
            _edge.nodeIndex2 = edgeIter->first;

            // if both nodes have high edge counts we filter out the edge:
            if (isEdgeFilterNode1) {
              const SVLocusNode& node2(locus.getNode(_edge.nodeIndex2));
              const bool         isEdgeFilterNode2(node2.size() > _graphNodeMaxEdgeCount);
              if (isEdgeFilterNode2) {
#ifdef DEBUG_EDGER
                log_os << "EDGER: jump filtering @ hc: " << _headCount << "\n";
#endif
                isLastFiltered = true;
                continue;
              }
            }
            return;
          }
        }

        _edge.nodeIndex1++;
      }

      // It is not common, but possible, to reach the end of the above while loop without returning.
      // This can happen if the final edge(s) in the locus are all filtered.
#ifdef DEBUG_EDGER
      log_os << "EDGER: reached end of locus without finding first edge: " << _edge.locusIndex << "\n";
#endif
      assert(isLastFiltered);
      assert(_headCount >= _beginCount);
    }
    _edge.locusIndex++;
    _edge.nodeIndex1 = 0;
    _edge.nodeIndex2 = 0;
  }

  assert(false && "jumpToFirstEdge: invalid state");
}

void EdgeRetrieverBin::advanceEdge()
{
  typedef SVLocusEdgesType::const_iterator edgeiter_t;

  const bool     isFilterNodes(_graphNodeMaxEdgeCount > 0);
  const unsigned setSize(_set.size());

  if (0 != _headCount) _edge.nodeIndex2++;

  bool isLastFiltered(false);

  while (true) {
    if (isLastFiltered && (_edge.locusIndex == setSize)) {
      _headCount       = (_endCount + 1);
      _edge.locusIndex = 0;
      _edge.nodeIndex1 = 0;
      _edge.nodeIndex2 = 0;
      return;
    }
    assert(_edge.locusIndex < setSize);

    const SVLocus& locus(_set.getLocus(_edge.locusIndex));
    const unsigned locusSize(locus.size());

    while (_edge.nodeIndex1 < locusSize) {
      const SVLocusNode&       node1(locus.getNode(_edge.nodeIndex1));
      const bool               isEdgeFilterNode1(isFilterNodes && (node1.size() > _graphNodeMaxEdgeCount));
      const SVLocusEdgeManager node1Manager(node1.getEdgeManager());
      edgeiter_t               edgeIter(node1Manager.getMap().lower_bound(_edge.nodeIndex2));
      const edgeiter_t         edgeIterEnd(node1Manager.getMap().cend());

      for (; edgeIter != edgeIterEnd; ++edgeIter) {
        unsigned   edgeCount(edgeIter->second.getCount());
        const bool isSelfEdge(edgeIter->first == _edge.nodeIndex1);
        if (!isSelfEdge) edgeCount += locus.getEdge(edgeIter->first, _edge.nodeIndex1).getCount();
        _headCount += edgeCount;
        _edge.nodeIndex2 = edgeIter->first;

        // if both nodes have high edge counts we filter out the edge:
        if (isEdgeFilterNode1) {
          const SVLocusNode& node2(locus.getNode(_edge.nodeIndex2));
          const bool         isEdgeFilterNode2(node2.size() > _graphNodeMaxEdgeCount);
          if (isEdgeFilterNode2) {
#ifdef DEBUG_EDGER
            log_os << "EDGER: advance filtering @ hc: " << _headCount << "\n";
#endif
            isLastFiltered = true;
            continue;
          }
        }

        return;
      }
      ++_edge.nodeIndex1;
      _edge.nodeIndex2 = _edge.nodeIndex1;
    }
    ++_edge.locusIndex;
    _edge.nodeIndex1 = 0;
    _edge.nodeIndex2 = 0;
  }

  assert(false && "advanceEdge: invalid state");
}

bool EdgeRetrieverBin::next()
{
#ifdef DEBUG_EDGER
  log_os << "EDGER: start next hc: " << _headCount << "\n";
#endif

  if (_headCount >= _endCount) return false;

  // first catch headCount up to the begin edge if required:
  if (_headCount < _beginCount) {
    jumpToFirstEdge();
#ifdef DEBUG_EDGER
    log_os << "EDGER: jumped hc: " << _headCount << " " << _edge << "\n";
#endif
  } else {
    advanceEdge();
#ifdef DEBUG_EDGER
    log_os << "EDGER: advanced hc: " << _headCount << " " << _edge << "\n";
#endif
  }

  assert(_headCount >= _beginCount);
  return (_headCount <= _endCount);
}
