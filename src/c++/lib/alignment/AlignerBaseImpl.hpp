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

#ifdef DEBUG_ALN_MATRIX
#include <boost/io/ios_state.hpp>

#include <iomanip>
#include <iostream>

template <typename ScoreType>
template <typename SymIter, typename MatrixType, typename ScoreValType>
void AlignerBase<ScoreType>::dumpSingleRefTable(
    const SymIter                                 refBegin,
    const SymIter                                 refEnd,
    const size_t                                  querySize,
    const MatrixType&                             ptrMatrix,
    const std::vector<std::vector<ScoreValType>>& storeScores,
    const char                                    refSym,
    const AlignState::index_t                     sIndex,
    unsigned&                                     storeIndex,
    std::ostream&                                 os) const
{
  boost::io::ios_all_saver guard(os);

  auto printVal = [](const ScoreType& val, const char fromSym, std::ostream& pos) {
    if (val < -900) {
      pos << " XX";
    } else {
      pos << std::setfill(' ') << std::setw(3) << val;
    }
    pos << fromSym;
  };

  auto printQueryRow = [&](const unsigned qrefIndex) {
    for (unsigned queryIndex(0); queryIndex <= querySize; ++queryIndex) {
      const auto& val(storeScores[storeIndex][queryIndex].getScore(sIndex));
      const char  fromSym(AlignState::symbol(ptrMatrix.val(queryIndex, qrefIndex).getStatePtr(sIndex)));
      printVal(val, fromSym, os);
    }
    os << "\n";
  };

  os << "# - ";
  printQueryRow(0);
  unsigned refIndex(0);
  for (SymIter refIter(refBegin); refIter != refEnd; ++refIter, ++refIndex) {
    os << refSym << " " << *refIter << " ";
    storeIndex++;
    printQueryRow(refIndex + 1);
  }
}
#endif
