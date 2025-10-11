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
