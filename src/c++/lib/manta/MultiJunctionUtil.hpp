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
