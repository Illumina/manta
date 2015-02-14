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
#include "svgraph/EdgeInfo.hh"

#include <cstdint>


/// aggregate statistics over a group of GSV edges
struct GSCEdgeGroupStats
{
    void
    merge(const GSCEdgeGroupStats& rhs)
    {
        totalTime += rhs.totalTime;
        assemblyTime += rhs.assemblyTime;
        scoringTime += rhs.scoringTime;
        total += rhs.total;
    }

    template<class Archive>
    void serialize(Archive& ar, const unsigned /* version */)
    {
        ar& totalTime& assemblyTime& scoringTime& total;
    }

    double totalTime = 0;
    double assemblyTime = 0;
    double scoringTime = 0;
    uint64_t total = 0;
};


struct GSCEdgeStats
{
    void
    merge(const GSCEdgeStats& rhs)
    {
        local.merge(rhs.local);
        remote.merge(rhs.remote);
    }

    void
    update(
        const EdgeInfo& edge,
        const EdgeRuntimeTracker& edgeTracker)
    {
        GSCEdgeGroupStats& gStats(getStatsGroup(edge));
        gStats.total += 1;
        gStats.totalTime += edgeTracker.getLastEdgeTime();
        gStats.assemblyTime += edgeTracker.assmTime.getSeconds();
        gStats.scoringTime += edgeTracker.scoreTime.getSeconds();
    }

    template<class Archive>
    void serialize(Archive& ar, const unsigned /* version */)
    {
        ar& local& remote;
    }

private:
    GSCEdgeGroupStats&
    getStatsGroup(const EdgeInfo& edge)
    {
        return (edge.isSelfEdge() ? local : remote);
    }

public:
    GSCEdgeGroupStats local;
    GSCEdgeGroupStats remote;
};



