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


/// parameters specific to SVLocusSet:
///
struct SVLocusSetOptions
{
    explicit
    SVLocusSetOptions(
        const unsigned initObservationWeight = 1) :
        observationWeight(initObservationWeight),
        minMergeEdgeObservations(3),
        maxSearchCount(500),
        maxSearchDensity(0.5)
    {}

    unsigned
    getMinMergeEdgeCount() const
    {
        return (observationWeight*minMergeEdgeObservations);
    }

    template<class Archive>
    void serialize(Archive& ar, const unsigned /* version */)
    {
        ar& observationWeight;
        ar& minMergeEdgeObservations;
        ar& maxSearchCount;
        ar& maxSearchDensity;
    }

    unsigned observationWeight; ///< used to translate graph edges counts to observations
    unsigned minMergeEdgeObservations; ///< to reduce noise in the graph, we only merge once shared edges reach this number of observations
    unsigned maxSearchCount; ///< the search for intersecting regions in the graph stops once this number is reached
    float maxSearchDensity; ///< the search for intersecting regions in the graph stops once this many regions/base are found
};

