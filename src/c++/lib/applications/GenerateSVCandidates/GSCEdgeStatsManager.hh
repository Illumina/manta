// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Manta
// Copyright (c) 2013-2015 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/sequencing/licenses/>
//

///
/// \author Chris Saunders
///

#pragma once

#include "EdgeRuntimeTracker.hh"
#include "appstats/GSCEdgeStats.hh"
#include "blt_util/time_util.hh"
#include "svgraph/EdgeInfo.hh"

#include "boost/utility.hpp"

#include <iosfwd>
#include <string>


/// handles all messy real world interaction for the stats module,
/// stats module itself just accumulates data and
///
struct GSCEdgeStatsManager : private boost::noncopyable
{
    GSCEdgeStatsManager(
        const std::string& outputFile);

    ~GSCEdgeStatsManager();

    void
    updateEdgeCandidates(
        const EdgeInfo& edge,
        const unsigned candCount,
        const SVFinderStats& finderStats)
    {
        if (_osPtr == nullptr) return;

        GSCEdgeGroupStats& gStats(getStatsGroup(edge));
        gStats.totalInputEdgeCount++;
        gStats.totalCandidateCount+=candCount;
        gStats.candidatesPerEdge.increment(candCount);
        gStats.finderStats.merge(finderStats);
    }

    void
    updateMJFilter(
        const EdgeInfo& edge,
        const unsigned mjComplexCount,
        const unsigned mjSpanningFilterCount)
    {
        if (_osPtr == nullptr) return;

        GSCEdgeGroupStats& gStats(getStatsGroup(edge));
        gStats.totalComplexCandidate += mjComplexCount;
        gStats.totalSpanningCandidateFilter += mjSpanningFilterCount;
    }

    void
    updateJunctionCandidates(
        const EdgeInfo& edge,
        const unsigned junctionCount,
        const bool isComplex)
    {
        if (_osPtr == nullptr) return;

        GSCEdgeGroupStats& gStats(getStatsGroup(edge));
        gStats.totalJunctionCount+=junctionCount;
        if (isComplex) gStats.totalComplexJunctionCount+=junctionCount;
        gStats.breaksPerJunction.increment(junctionCount);
    }

    void
    updateAssemblyCount(
        const EdgeInfo& edge,
        const unsigned assemblyCount,
        const bool isSpanning,
        const bool isOverlapSkip = false)
    {
        if (_osPtr == nullptr) return;

        GSCEdgeGroupStats& gStats(getStatsGroup(edge));
        gStats.totalAssemblyCandidates += assemblyCount;
        if (isSpanning) gStats.totalSpanningAssemblyCandidates += assemblyCount;
        if (isOverlapSkip)
        {
            gStats.totalJunctionAssemblyOverlapSkips++;
        }
        else
        {
            gStats.assemblyCandidatesPerJunction.increment(assemblyCount);
        }
    }


    void
    updateScoredEdgeTime(
        const EdgeInfo& edge,
        const EdgeRuntimeTracker& edgeTracker)
    {
        if (_osPtr == nullptr) return;

        GSCEdgeGroupStats& gStats(getStatsGroup(edge));
        gStats.totalTime.merge(edgeTracker.getLastEdgeTime());
        gStats.candTime.merge(edgeTracker.candTime.getTimes());
        gStats.assemblyTime.merge(edgeTracker.assmTime.getTimes());
        gStats.scoringTime.merge(edgeTracker.scoreTime.getTimes());
    }

private:
    GSCEdgeGroupStats&
    getStatsGroup(
        const EdgeInfo& edge)
    {
        return (edge.isSelfEdge() ? edgeStats.edgeData.selfEdges : edgeStats.edgeData.remoteEdges);
    }

    std::ostream* _osPtr;
    TimeTracker lifeTime;
    GSCEdgeStats edgeStats;
};
