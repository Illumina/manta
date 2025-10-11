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

#pragma once

template <typename ScoreType>
struct AlignmentScores {
  AlignmentScores(
      ScoreType initMatch,
      ScoreType initMismatch,
      ScoreType initOpen,
      ScoreType initExtend,
      ScoreType initOffEdge,
      bool      initIsAllowEdgeInsertion = false)
    : match(initMatch),
      mismatch(initMismatch),
      open(initOpen),
      extend(initExtend),
      offEdge(initOffEdge),
      isAllowEdgeInsertion(initIsAllowEdgeInsertion)
  {
  }

  /// Match score
  const ScoreType match;

  /// Mismatch score (should be negative)
  const ScoreType mismatch;

  /// Gap open, gap of length N is scored (open + N * extend) (should be negative)
  const ScoreType open;

  /// Gap extend, gap of length N is scored (open + N * extend) (should be negative or zero)
  const ScoreType extend;

  /// Score applied when query goes off the end of an edge (should be negative)
  const ScoreType offEdge;

  /// Are insertions allowed directly on the leading and trailing edges?
  const bool isAllowEdgeInsertion;
};
