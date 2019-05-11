//
// Manta - Structural Variant and Indel Caller
// Copyright (c) 2013-2019 Illumina, Inc.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//

/// \file
/// \brief Utility functions for mocking SVLocus elements during unit testing
/// \author Trevor Ramsay
///

#pragma once

#include "svgraph/SVLocus.hpp"

/// \brief Add a pair of nodes to a SVLocus object.
/// \param count The evidence count applied to the edge from node1 to node2
void locusAddPair(
    SVLocus&       locus,
    const int32_t  tid1,
    const int32_t  beginPos1,
    const int32_t  endPos1,
    const int32_t  tid2,
    const int32_t  beginPos2,
    const int32_t  endPos2,
    const bool     bothLocal = false,
    const unsigned count     = 1);
