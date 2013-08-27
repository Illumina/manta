// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Manta
// Copyright (c) 2013 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/downloads/sequencing/licenses/>.
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
        ar& minMergeEdgeCount& maxSearchCount& maxSearchDensity;
    }

    unsigned minMergeEdgeCount; ///< to reduce noise in the graph, we only merge once shared edges (including self-edges) reach this count
    unsigned maxSearchCount;
    float maxSearchDensity;
};

