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
