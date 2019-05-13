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

#include <cassert>

//#define DEBUG_ALN

#ifdef DEBUG_ALN
#include <iostream>
#include "blt_util/log.hpp"
#endif

template <typename ScoreType>
template <typename SymIter>
void GlobalLargeIndelAligner<ScoreType>::align(
    const SymIter               queryBegin,
    const SymIter               queryEnd,
    const SymIter               refBegin,
    const SymIter               refEnd,
    AlignmentResult<ScoreType>& result) const
{
  result.clear();

  const AlignmentScores<ScoreType>& scores(this->getScores());

  const size_t querySize(std::distance(queryBegin, queryEnd));
  const size_t refSize(std::distance(refBegin, refEnd));

  assert(0 != querySize);
  assert(0 != refSize);

  _score1.resize(querySize + 1);
  _score2.resize(querySize + 1);
  _ptrMat.resize(querySize + 1, refSize + 1);

  ScoreVec* thisSV(&_score1);
  ScoreVec* prevSV(&_score2);

  static const ScoreType badVal(-10000);

  // global alignment of query
  //
  // disallow start from the delete state, control start from insert state with flag
  //
  // query can 'fall-off' the end of a short reference, in which case it will
  // be soft-clipped and each base off the end will be scored as offEdge
  //
  for (unsigned queryIndex(0); queryIndex <= querySize; queryIndex++) {
    PtrVal&   headPtr(_ptrMat.val(queryIndex, 0));
    ScoreVal& val((*thisSV)[queryIndex]);
    headPtr.match = AlignState::MATCH;
    val.match     = queryIndex * scores.offEdge;
    headPtr.del   = AlignState::MATCH;
    val.del       = badVal;
    if (not scores.isAllowEdgeInsertion) {
      headPtr.ins = AlignState::MATCH;
      val.ins     = badVal;
    } else {
      headPtr.ins = AlignState::INSERT;
      val.ins     = scores.open + (queryIndex * scores.extend);
    }
    headPtr.jumpDel = AlignState::MATCH;
    val.jumpDel     = badVal;
    headPtr.jumpIns = AlignState::MATCH;
    val.jumpIns     = badVal;
  }

#ifdef DEBUG_ALN_MATRIX
  // store full matrix of scores to print out later, don't turn this debug option on for large references!
  std::vector<ScoreVec> storeScores;

  storeScores.push_back(*thisSV);
#endif

  BackTrace<ScoreType> btrace;

  {
    unsigned refIndex(0);
    for (SymIter refIter(refBegin); refIter != refEnd; ++refIter, ++refIndex) {
      std::swap(thisSV, prevSV);

      {
        // disallow start from the insert or delete states:
        PtrVal&   headPtr(_ptrMat.val(0, refIndex + 1));
        ScoreVal& val((*thisSV)[0]);
        headPtr.match   = AlignState::MATCH;
        val.match       = 0;
        headPtr.del     = AlignState::MATCH;
        val.del         = badVal;
        headPtr.ins     = AlignState::MATCH;
        val.ins         = badVal;
        headPtr.jumpDel = AlignState::MATCH;
        val.jumpDel     = badVal;
        headPtr.jumpIns = AlignState::MATCH;
        val.jumpIns     = badVal;
      }

      unsigned queryIndex(0);
      for (SymIter queryIter(queryBegin); queryIter != queryEnd; ++queryIter, ++queryIndex) {
        // update match
        ScoreVal& headScore((*thisSV)[queryIndex + 1]);
        PtrVal&   headPtr(_ptrMat.val(queryIndex + 1, refIndex + 1));
        {
          const ScoreVal& sval((*prevSV)[queryIndex]);
          headPtr.match =
              this->max5(headScore.match, sval.match, sval.del, sval.ins, sval.jumpDel, sval.jumpIns);

          headScore.match += ((*queryIter == *refIter) ? scores.match : scores.mismatch);
        }

        // update delete
        {
          const ScoreVal& sval((*prevSV)[queryIndex + 1]);
          headPtr.del =
              this->max5(headScore.del, sval.match + scores.open, sval.del, sval.ins, badVal, sval.jumpIns);

          headScore.del += scores.extend;
          if (0 == queryIndex) headScore.del = badVal;
        }

        // update insert
        {
          const ScoreVal& sval((*thisSV)[queryIndex]);
          headPtr.ins = this->max5(headScore.ins, sval.match + scores.open, badVal, sval.ins, badVal, badVal);

          headScore.ins += scores.extend;
          if (0 == queryIndex) headScore.ins = badVal;
        }

        // update jumpDel
        {
          // you can switch between long ins and delete but only
          // by paying the full open penalty again
          //
          // the switch from short ins to jump del is meant to simulate
          // a free transition from jump del to short ins, but makes
          // the cigar I->D order come out the same as other aligners in this library
          const ScoreVal& sval((*prevSV)[queryIndex + 1]);
          headPtr.jumpDel = this->max5(
              headScore.jumpDel,
              sval.match + _largeIndelScore,
              badVal,
              sval.ins + _largeIndelScore - scores.open,
              sval.jumpDel,
              sval.jumpIns + _largeIndelScore);

          if (0 == queryIndex) headScore.jumpDel = badVal;
        }

        // update jumpIns
        {
          const ScoreVal& sval((*thisSV)[queryIndex]);
          headPtr.jumpIns = this->max5(
              headScore.jumpIns, sval.match + _largeIndelScore, badVal, badVal, badVal, sval.jumpIns);

          if (0 == queryIndex) headScore.jumpIns = badVal;
        }

#ifdef DEBUG_ALN
        log_os << "queryIdx refIdx: " << queryIndex + 1 << " " << refIndex + 1 << "\n";
        log_os << headScore.match << ":" << headScore.del << ":" << headScore.ins << ":" << headScore.jumpDel
               << ":" << headScore.jumpIns << "/" << static_cast<int>(headPtr.match)
               << static_cast<int>(headPtr.del) << static_cast<int>(headPtr.ins)
               << static_cast<int>(headPtr.jumpDel) << static_cast<int>(headPtr.jumpIns) << "\n";
#endif
      }
#ifdef DEBUG_ALN
      log_os << "\n";
#endif

#ifdef DEBUG_ALN_MATRIX
      storeScores.push_back(*thisSV);
#endif

      // get backtrace info:
      {
        const ScoreVal& sval((*thisSV)[querySize]);
        updateBacktrace(sval.match, refIndex + 1, querySize, btrace);
      }
    }
  }

  // optionally allow for trailing insertion
  if (scores.isAllowEdgeInsertion) {
    const ScoreVal& sval((*thisSV)[querySize]);
    updateBacktrace(sval.ins, refSize, querySize, btrace, AlignState::INSERT);
  }

  // in the backtrace start search, also allow for the case where the query falls-off the end of the
  // reference:
  for (unsigned queryIndex(0); queryIndex <= querySize; queryIndex++) {
    const ScoreVal& sval((*thisSV)[queryIndex]);
    const ScoreType thisMax(sval.match + (querySize - queryIndex) * scores.offEdge);
    updateBacktrace(thisMax, refSize, queryIndex, btrace);
  }

#ifdef DEBUG_ALN_MATRIX
  std::vector<AlignState::index_t> dumpStates{
      AlignState::MATCH, AlignState::DELETE, AlignState::INSERT, AlignState::JUMP, AlignState::JUMPINS};
  this->dumpTables(queryBegin, queryEnd, refBegin, refEnd, querySize, _ptrMat, dumpStates, storeScores);
#endif

  this->backTraceAlignment(
      queryBegin, queryEnd, refBegin, refEnd, querySize, refSize, _ptrMat, btrace, result);
}
