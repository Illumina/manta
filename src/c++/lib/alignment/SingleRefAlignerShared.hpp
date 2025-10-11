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
/// \brief Common boilerplate for single-reference sequence aligners

#pragma once

#include "AlignerBase.hpp"
#include "AlignerUtil.hpp"
#include "Alignment.hpp"

#include "blt_util/basic_matrix.hpp"

#include <iosfwd>

template <typename ScoreType>
struct AlignmentResult {
  AlignmentResult() { clear(); }

  void clear()
  {
    score    = 0;
    isJumped = false;
    align.clear();
  }

  ScoreType score;
  bool      isJumped;  ///< whether alignment path includes jump state(s) while backtracking
  Alignment align;
};

template <typename ScoreType>
std::ostream& operator<<(std::ostream& os, AlignmentResult<ScoreType>& alignment);

template <typename ScoreType>
struct SingleRefAlignerBase : public AlignerBase<ScoreType> {
  SingleRefAlignerBase(const AlignmentScores<ScoreType>& scores) : AlignerBase<ScoreType>(scores) {}

protected:
  /// returns alignment path of query to reference
  template <typename SymIter, typename MatrixType>
  void backTraceAlignment(
      const SymIter               queryBegin,
      const SymIter               queryEnd,
      const SymIter               refBegin,
      const SymIter               refEnd,
      const size_t                querySize,
      const size_t                refSize,
      const MatrixType&           ptrMatrix,
      const BackTrace<ScoreType>& btraceInput,
      AlignmentResult<ScoreType>& result) const;

#ifdef DEBUG_ALN_MATRIX
  /// write out matrix of scores and back-trace pointers for debug:
  template <typename SymIter, typename MatrixType, typename ScoreValType>
  void dumpTables(
      const SymIter                                 queryBegin,
      const SymIter                                 queryEnd,
      const SymIter                                 refBegin,
      const SymIter                                 refEnd,
      const size_t                                  querySize,
      const MatrixType&                             ptrMatrix,
      const std::vector<AlignState::index_t>&       dumpStates,
      const std::vector<std::vector<ScoreValType>>& storeScores) const;
#endif
};

#include "SingleRefAlignerSharedImpl.hpp"
