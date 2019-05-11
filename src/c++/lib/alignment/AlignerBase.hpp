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

#include "alignment/Alignment.hpp"
#include "alignment/AlignmentScores.hpp"

// #define DEBUG_ALN // Standard debug output
// #define DEBUG_ALN_MATRIX // Dump full edit-matrix tables to stderr. Does not scale to non-trivial ref/query size!

#ifdef DEBUG_ALN_MATRIX
#include <iosfwd>
#endif

/// shared methods for all aligners
///
template <typename ScoreType>
struct AlignerBase {
  AlignerBase(const AlignmentScores<ScoreType>& scores) : _scores(scores) {}

  /// read-only access to the aligner's scores:
  const AlignmentScores<ScoreType>& getScores() const { return _scores; }

protected:
  static uint8_t max3(ScoreType& max, const ScoreType v0, const ScoreType v1, const ScoreType v2)
  {
    max         = v0;
    uint8_t ptr = 0;
    if (v1 > v0) {
      max = v1;
      ptr = 1;
    }
    if (v2 > max) {
      max = v2;
      ptr = 2;
    }
    return ptr;
  }

#ifdef DEBUG_ALN_MATRIX
  /// write out subset of matrix of scores back-trace pointers for debug,for
  /// one reference
  template <typename SymIter, typename MatrixType, typename ScoreValType>
  void dumpSingleRefTable(
      const SymIter                                 refBegin,
      const SymIter                                 refEnd,
      const size_t                                  querySize,
      const MatrixType&                             ptrMatrix,
      const std::vector<std::vector<ScoreValType>>& storeScores,
      const char                                    refSym,
      const AlignState::index_t                     sIndex,
      unsigned&                                     storeIndex,
      std::ostream&                                 os) const;
#endif

  const AlignmentScores<ScoreType> _scores;
};

#include "alignment/AlignerBaseImpl.hpp"
