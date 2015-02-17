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
        totalEdgeCount += rhs.totalEdgeCount;
    }

    template<class Archive>
    void serialize(Archive& ar, const unsigned /* version */)
    {
        ar& BOOST_SERIALIZATION_NVP(totalTime)& BOOST_SERIALIZATION_NVP(assemblyTime)& BOOST_SERIALIZATION_NVP(scoringTime)& BOOST_SERIALIZATION_NVP(totalEdgeCount);
    }

    double totalTime = 0;
    double assemblyTime = 0;
    double scoringTime = 0;
    uint64_t totalEdgeCount = 0;
};

BOOST_CLASS_IMPLEMENTATION(GSCEdgeGroupStats, boost::serialization::object_serializable)


struct GSCEdgeStats
{
    GSCEdgeStats() {}

    void
    merge(const GSCEdgeStats& rhs)
    {
        selfEdges.merge(rhs.selfEdges);
        remoteEdges.merge(rhs.remoteEdges);
    }

    template<class Archive>
    void serialize(Archive& ar, const unsigned /* version */)
    {
        ar& BOOST_SERIALIZATION_NVP(selfEdges) & BOOST_SERIALIZATION_NVP(remoteEdges);
    }

    GSCEdgeGroupStats selfEdges;
    GSCEdgeGroupStats remoteEdges;
};

BOOST_CLASS_IMPLEMENTATION(GSCEdgeStats, boost::serialization::object_serializable)


