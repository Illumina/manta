// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Manta
// Copyright (c) 2013-2014 Illumina, Inc.
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
#include "blt_util/time_util.hh"
#include "svgraph/GSCEdgeStats.hh"
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
        const unsigned candCount)
    {
        if (_osPtr == nullptr) return;

        GSCEdgeGroupStats& gStats(getStatsGroup(edge));
        gStats.totalInputEdgeCount++;
        gStats.totalCandidateCount+=candCount;
        gStats.candidatesPerEdge.increment(candCount);
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
        const unsigned assemblyCount = 0,
        const bool isSpanning = false)
    {
        if (_osPtr == nullptr) return;

        GSCEdgeGroupStats& gStats(getStatsGroup(edge));
        gStats.totalAssemblyCandidates += assemblyCount;
        if (isSpanning) gStats.totalSpanningAssemblyCandidates += assemblyCount;
        gStats.assemblyCandidatesPerJunction.increment(assemblyCount);
    }


    void
    updateScoredEdgeTime(
        const EdgeInfo& edge,
        const EdgeRuntimeTracker& edgeTracker)
    {
        if (_osPtr == nullptr) return;

        GSCEdgeGroupStats& gStats(getStatsGroup(edge));
        gStats.totalTime += edgeTracker.getLastEdgeTime();
        gStats.candTime += edgeTracker.candTime.getSeconds();
        gStats.assemblyTime += edgeTracker.assmTime.getSeconds();
        gStats.scoringTime += edgeTracker.scoreTime.getSeconds();
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
