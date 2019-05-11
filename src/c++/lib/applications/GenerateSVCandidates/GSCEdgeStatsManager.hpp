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

#pragma once

#include "EdgeRuntimeTracker.hpp"
#include "appstats/GSCEdgeStats.hpp"
#include "blt_util/time_util.hpp"
#include "svgraph/EdgeInfo.hpp"

#include "boost/utility.hpp"

#include <iosfwd>
#include <string>

/// This object processes incoming edge data to accumulate stats in the edge stats module (GSCEdgeStats)
///
struct GSCEdgeStatsManager : private boost::noncopyable {
  explicit GSCEdgeStatsManager() { lifeTime.resume(); }

  GSCEdgeStats returnStats()
  {
    edgeStats.edgeData.lifeTime = lifeTime.getTimes();
    return edgeStats;
  }

  void updateEdgeCandidates(const EdgeInfo& edge, const unsigned candCount, const SVFinderStats& finderStats)
  {
    GSCEdgeGroupStats& gStats(getStatsGroup(edge));
    gStats.totalInputEdgeCount++;
    gStats.totalCandidateCount += candCount;
    gStats.candidatesPerEdge.increment(candCount);
    gStats.finderStats.merge(finderStats);
  }

  void updateMJFilter(
      const EdgeInfo& edge, const unsigned mjComplexCount, const unsigned mjSpanningFilterCount)
  {
    GSCEdgeGroupStats& gStats(getStatsGroup(edge));
    gStats.totalComplexCandidate += mjComplexCount;
    gStats.totalSpanningCandidateFilter += mjSpanningFilterCount;
  }

  void updateJunctionCandidateCounts(const EdgeInfo& edge, const unsigned junctionCount, const bool isComplex)
  {
    GSCEdgeGroupStats& gStats(getStatsGroup(edge));
    gStats.totalJunctionCount += junctionCount;
    if (isComplex) gStats.totalComplexJunctionCount += junctionCount;
    gStats.breaksPerJunction.increment(junctionCount);
  }

  void updateAssemblyCount(
      const EdgeInfo& edge,
      const unsigned  assemblyCount,
      const bool      isSpanning,
      const bool      isOverlapSkip = false)
  {
    GSCEdgeGroupStats& gStats(getStatsGroup(edge));
    gStats.totalAssemblyCandidates += assemblyCount;
    if (isSpanning) gStats.totalSpanningAssemblyCandidates += assemblyCount;
    if (isOverlapSkip) {
      gStats.totalJunctionAssemblyOverlapSkips++;
    } else {
      gStats.assemblyCandidatesPerJunction.increment(assemblyCount);
    }
  }

  void updateScoredEdgeTime(const EdgeInfo& edge, const EdgeRuntimeTracker& edgeTracker)
  {
    GSCEdgeGroupStats& gStats(getStatsGroup(edge));
    gStats.totalTime.merge(edgeTracker.getLastEdgeTime());
    gStats.candTime.merge(edgeTracker.candidacyTime.getTimes());
    gStats.assemblyTime.merge(edgeTracker.assemblyTime.getTimes());
    gStats.scoringTime.merge(edgeTracker.scoreTime.getTimes());
  }

private:
  GSCEdgeGroupStats& getStatsGroup(const EdgeInfo& edge)
  {
    return (edge.isSelfEdge() ? edgeStats.edgeData.selfEdges : edgeStats.edgeData.remoteEdges);
  }

  /// Lifetime of thread, whether it's running or not
  TimeTracker  lifeTime;
  GSCEdgeStats edgeStats;
};
