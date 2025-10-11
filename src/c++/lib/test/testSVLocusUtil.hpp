//
// Manta - Structural Variant and Indel Caller
// Copyright (c) 2013-2025 Illumina, Inc.
//
// This program is licensed under the terms of the Polyform strict license
//
// ***As far as the law allows, the software comes as is, without
// any warranty or condition, and the licensor will not be liable
// to you for any damages arising out of these terms or the use
// or nature of the software, under any kind of legal claim.***
//
// You should have received a copy of the PolyForm Strict License 1.0.0
// along with this program.  If not, see <https://polyformproject.org/licenses/strict/1.0.0>.
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
