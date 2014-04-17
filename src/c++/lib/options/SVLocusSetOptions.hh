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


/// parameters specific to SVLocusSet:
///
struct SVLocusSetOptions
{
    SVLocusSetOptions(
        const unsigned initMinMergeEdgeCount = 2,
        const unsigned initMaxSearchCount = 500,
        const float initMaxSearchDensity = 0.5) :
        minMergeEdgeCount(initMinMergeEdgeCount),
        maxSearchCount(initMaxSearchCount),
        maxSearchDensity(initMaxSearchDensity)
    {}

    template<class Archive>
    void serialize(Archive& ar, const unsigned /* version */)
    {
        ar& minMergeEdgeCount;
        ar& maxSearchCount;
        ar& maxSearchDensity;
    }

    unsigned minMergeEdgeCount; ///< to reduce noise in the graph, we only merge once shared edges reach this count
    unsigned maxSearchCount; ///< the search for intersecting regions in the graph stops once this number is reached
    float maxSearchDensity; ///< the search for intersecting regions in the graph stops once this many regions/base are found
};

