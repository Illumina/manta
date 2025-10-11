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

#include <iostream>

#include <cassert>

//#define DEBUG_ALN

#if defined(DEBUG_ALN) || defined(DEBUG_ALN_MATRIX)
#include "blt_util/log.hpp"
#endif

#ifdef DEBUG_ALN_MATRIX
template <typename ScoreType>
template <typename SymIter, typename MatrixType, typename ScoreValType>
void SingleRefAlignerBase<ScoreType>::dumpTables(
    const SymIter                                 queryBegin,
    const SymIter                                 queryEnd,
    const SymIter                                 refBegin,
    const SymIter                                 refEnd,
    const size_t                                  querySize,
    const MatrixType&                             ptrMatrix,
    const std::vector<AlignState::index_t>&       dumpStates,
    const std::vector<std::vector<ScoreValType>>& storeScores) const
{
  const unsigned stateSize(dumpStates.size());
  for (unsigned stateIndex(0); stateIndex < stateSize; ++stateIndex) {
    const AlignState::index_t sIndex(dumpStates[stateIndex]);
    log_os << "******** Dumping matrix for state: " << AlignState::label(sIndex) << " ********\n";
    {
      unsigned queryIndex(0);
      log_os << "REF    -";
      for (SymIter queryIter(queryBegin); queryIter != queryEnd; ++queryIter, ++queryIndex) {
        log_os << "   " << *queryIter;
      }
      log_os << "\n";
    }

    // dump state before refIndex 0
    unsigned storeIndex(0);
    this->dumpSingleRefTable(
        refBegin, refEnd, querySize, ptrMatrix, storeScores, '1', sIndex, storeIndex, log_os);
  }
}
#endif

template <typename ScoreType>
std::ostream& operator<<(std::ostream& os, AlignmentResult<ScoreType>& alignment)
{
  os << "AlignerResult: score: " << alignment.score << "\n"
     << "\talignment: " << alignment.align << "\n";
  return os;
}

// traceback:
template <typename ScoreType>
template <typename SymIter, typename MatrixType>
void SingleRefAlignerBase<ScoreType>::backTraceAlignment(
    const SymIter               queryBegin,
    const SymIter               queryEnd,
    const SymIter               refBegin,
    const SymIter               refEnd,
    const size_t                querySize,
    const size_t                refSize,
    const MatrixType&           ptrMatrix,
    const BackTrace<ScoreType>& btraceInput,
    AlignmentResult<ScoreType>& result) const
{
  BackTrace<ScoreType> btrace(btraceInput);

  assert(btrace.isInit);
  assert(btrace.refBegin <= refSize);
  assert(btrace.queryBegin <= querySize);

#ifdef DEBUG_ALN
  log_os << "qSize: " << querySize << " refSize: " << refSize << "\n";
  log_os << "bt-start max: " << btrace.max << " refBegin: " << btrace.refBegin
         << " qBegin: " << btrace.queryBegin << " state: " << AlignState::label(btrace.state) << "\n";
#endif

  result.score = btrace.max;

  // traceback:
  ALIGNPATH::path_t&      apath(result.align.apath);
  ALIGNPATH::path_segment ps;

  // add any trailing soft-clip if we go off the end of the reference:
  if (btrace.queryBegin < querySize) {
    ps.type   = ALIGNPATH::SOFT_CLIP;
    ps.length = (querySize - btrace.queryBegin);
  }

  while (true) {
    const AlignState::index_t nextState(
        ptrMatrix.val(btrace.queryBegin, btrace.refBegin).getStatePtr(btrace.state));

#ifdef DEBUG_ALN
    log_os << "bt-iter queryIndex: " << btrace.queryBegin << " refIndex: " << btrace.refBegin
           << " state: " << AlignState::label(btrace.state) << " next: " << AlignState::label(nextState)
           << "\n";
#endif

    if (btrace.state == AlignState::MATCH) {
      if ((btrace.queryBegin < 1) or (btrace.refBegin < 1)) break;
      AlignerUtil::updatePath(apath, ps, ALIGNPATH::MATCH);
      btrace.queryBegin--;
      btrace.refBegin--;
    } else if ((btrace.state == AlignState::DELETE) || (btrace.state == AlignState::JUMP)) {
      if (btrace.refBegin < 1) break;
      AlignerUtil::updatePath(apath, ps, ALIGNPATH::DELETE);
      btrace.refBegin--;
    } else if ((btrace.state == AlignState::INSERT) || (btrace.state == AlignState::JUMPINS)) {
      if (btrace.queryBegin < 1) break;
      AlignerUtil::updatePath(apath, ps, ALIGNPATH::INSERT);
      btrace.queryBegin--;
    } else {
      assert(false && "Unknown align state");
    }

    // check if the alignment path includes JUMP or JUMPINS states
    if ((btrace.state == AlignState::JUMP) || (btrace.state == AlignState::JUMPINS)) {
      result.isJumped = true;
#ifdef DEBUG_ALN
      log_os << "isJumped is set true"
             << "\n";
#endif
    }

    btrace.state = nextState;
    ps.length++;
  }

  if (ps.type != ALIGNPATH::NONE) apath.push_back(ps);

  // soft-clip beginning of read if we fall off the end of the reference
  if (btrace.queryBegin != 0) {
    ps.type   = ALIGNPATH::SOFT_CLIP;
    ps.length = btrace.queryBegin;
    apath.push_back(ps);
  }

  result.align.beginPos = btrace.refBegin;
  std::reverse(apath.begin(), apath.end());

  // if true, output final cigars using seq match '=' and mismatch 'X' symbols:
  static const bool isOutputSeqMatch(true);

  if (isOutputSeqMatch) {
    apath_add_seqmatch(queryBegin, queryEnd, (refBegin + result.align.beginPos), refEnd, apath);
  }
}
