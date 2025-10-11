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
