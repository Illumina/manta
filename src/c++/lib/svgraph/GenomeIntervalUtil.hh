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

#include "svgraph/GenomeInterval.hh"

#include <vector>


/// given a collection of genome intervals, reduce down to the minimum non-overlapping set:
///
/// \returns a vector with size equal to the input vector, containing a mapping of the input
///         interval index to the output interval index
///
std::vector<unsigned>
intervalCompressor(
    std::vector<GenomeInterval>& intervals);
