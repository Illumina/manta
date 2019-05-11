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
