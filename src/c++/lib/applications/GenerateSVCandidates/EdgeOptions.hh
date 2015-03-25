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


/// options for SVLocusGraph edge iteration and noise edge filtration
struct LocusEdgeOptions
{
    unsigned locusIndex = 0; ///< if isLocusIndex, report this locus only
    bool isNodeIndex1 = false; ///< if true, generate candidates for all edges touching a specifc node in one locus. Assumes isLocusIndex is true
    unsigned nodeIndex1 = 0;
    bool isNodeIndex2 = false; ///< if true, generate candidates for only the edge from node1 to node2 in one locus. Assumes isLocusIndex & isNodeIndex1 are true
    unsigned nodeIndex2 = 0;
};


/// options for SVLocusGraph edge iteration and noise edge filtration
struct EdgeOptions
{
    unsigned binCount = 1; ///< divide all edges in the graph into binCount bins of approx equal complexity
    unsigned binIndex = 0; ///< out of binCount bins, iterate through the edges in this bin only

    bool isLocusIndex = false; ///< if true, generate candidates for a specific SVgraph locus only, and ignore binCount/binIndex
    LocusEdgeOptions locusOpt;

    unsigned graphNodeMaxEdgeCount = 10; ///< if both nodes of an edge have an edge count higher than this, then skip evaluation of this edge, set to 0 to turn this filtration off
};
