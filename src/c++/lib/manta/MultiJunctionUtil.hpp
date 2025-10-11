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
/// \author Chris Saunders
///

#pragma once

#include <vector>

#include "manta/SVCandidate.hpp"
#include "manta/SVMultiJunctionCandidate.hpp"

/// \brief Convert independent SV candidates into multi-junction event candidates
///
/// Given a set of un-associated single-junction SV candidates, analyze which
/// candidate junctions could potentially be treated as a single multi-junction
/// event (such as a reciprocal translocation)
///
/// Right now multi-junction events are limited to pairs of (spannning) SV candidates, where
/// the breakends of both junctions in the pair are proximal and meeting the expected orientation
/// pattern consistent with a reciprocal translocation.
///
/// Note that 'complex' (short assembly targets from self-edges) SVs will not be grouped into
//
void findMultiJunctionCandidates(
    const std::vector<SVCandidate>&        svs,
    const unsigned                         minCandidateSpanningCount,
    const bool                             isRNA,
    unsigned&                              mjComplexCount,
    unsigned&                              mjSpanningFilterCount,
    std::vector<SVMultiJunctionCandidate>& mjSVs);
