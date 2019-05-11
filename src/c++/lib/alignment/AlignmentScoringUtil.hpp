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

#include "alignment/AlignmentScores.hpp"
#include "blt_util/align_path.hpp"

/// Get the alignment score for a given set of alignment scoring weights and alignment path
///
/// This requires SEQ_MATCH style alignment path which indicates match vs. mismatch positions
///
template <typename ScoreType>
ScoreType getPathScore(
    const AlignmentScores<ScoreType>& scores,
    const ALIGNPATH::path_t&          apath,
    const bool                        isScoreOffEdge = false);

/// Get the maximum partial path alignment score for a given set of alignment scoring weights and alignment
/// path
///
/// This requires SEQ_MATCH style alignment path which indicates match vs. mismatch positions
///
template <typename ScoreType>
ScoreType getMaxPathScore(
    const AlignmentScores<ScoreType>& scores,
    const ALIGNPATH::path_t&          apath,
    unsigned&                         maxReadOffset,
    unsigned&                         maxRefOffset,
    const bool                        isScoreOffEdge = true);

#include "alignment/AlignmentScoringUtilImpl.hpp"
