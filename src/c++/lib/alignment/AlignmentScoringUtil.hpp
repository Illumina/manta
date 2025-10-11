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
