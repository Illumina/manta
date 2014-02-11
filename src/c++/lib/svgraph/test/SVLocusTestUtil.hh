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

#include "svgraph/SVLocus.hh"


inline
void
locusAddPair(
    SVLocus& locus,
    const int32_t tid1,
    const int32_t beginPos1,
    const int32_t endPos1,
    const int32_t tid2,
    const int32_t beginPos2,
    const int32_t endPos2,
    const bool bothLocal = false,
    const int count = 1)
{
    const NodeIndexType nodePtr1 = locus.addNode(GenomeInterval(tid1,beginPos1,endPos1));
    const NodeIndexType nodePtr2 = locus.addNode(GenomeInterval(tid2,beginPos2,endPos2));
    const int remoteCount(bothLocal ? count : 0);
    locus.linkNodes(nodePtr1,nodePtr2,count,remoteCount);
}

