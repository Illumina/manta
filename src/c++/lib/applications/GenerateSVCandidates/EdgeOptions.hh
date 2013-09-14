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


struct EdgeOptions
{
    EdgeOptions() :
        binCount(1),
        binIndex(0),
        isLocusIndex(false),
        locusIndex(0),
        graphNodeMaxEdgeCount(10)
    {}

    unsigned binCount;
    unsigned binIndex;

    bool isLocusIndex; ///< if true, generate candidates for a specific SVgraph locus only
    unsigned locusIndex;

    unsigned graphNodeMaxEdgeCount; ///< if both nodes of an edge have an edge count higher than this, then skip evaluation of this edge, set to 0 to turn this filtration off
};

