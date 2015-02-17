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

#include "boost/serialization/nvp.hpp"

#include <cstdint>


/// aggregate statistics over a group of GSV edges
struct GSCEdgeGroupStats
{
    GSCEdgeGroupStats() {}

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
        ar& BOOST_SERIALIZATION_NVP(totalTime)& BOOST_SERIALIZATION_NVP(assemblyTime)& BOOST_SERIALIZATION_NVP(scoringTime)& BOOST_SERIALIZATION_NVP(total);
    }

    double totalTime = 0;
    double assemblyTime = 0;
    double scoringTime = 0;
    uint64_t total = 0;
};


struct GSCEdgeStats
{
    GSCEdgeStats() {}

    void
    merge(const GSCEdgeStats& rhs)
    {
        local.merge(rhs.local);
        remote.merge(rhs.remote);
    }

    template<class Archive>
    void serialize(Archive& ar, const unsigned /* version */)
    {
        ar& BOOST_SERIALIZATION_NVP(local) & BOOST_SERIALIZATION_NVP(remote);
    }

    GSCEdgeGroupStats local;
    GSCEdgeGroupStats remote;
};
