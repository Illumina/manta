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


/// aggregate statistics over a group of GSV edges
struct GSVEdgeGroupStats
{
    void
    merge(const GSVEdgeGroupStats& rhs)
    {
        candidateTime += rhs.candidateTime;
        assemblyTime += rhs.assemblyTime;
        scoringTime += rhs.scoringTime;
        total += rhs.total;
    }

    template<class Archive>
    void serialize(Archive& ar, const unsigned /* version */)
    {
        ar& candidateTime& assemblyTime& scoringTime& total;
    }

    double candidateTime = 0;
    double assemblyTime = 0;
    double scoringTime = 0;
    uint64_t total = 0;
};


struct GSVEdgeStats
{
    void
    merge(const GSVEdgeStats& rhs)
    {
        local.merge(rhs.local);
        remote.merge(rhs.remote);
    }

    template<class Archive>
    void serialize(Archive& ar, const unsigned /* version */)
    {
        ar& local& remote;
    }

    GSVEdgeGroupStats local;
    GSVEdgeGroupStats remote;
};
