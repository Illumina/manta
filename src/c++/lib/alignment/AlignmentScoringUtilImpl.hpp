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

//#define DEBUG_PATHSCORE

#ifdef DEBUG_PATHSCORE
#include "blt_util/log.hpp"

#include <iostream>
#endif

template <typename ScoreType>
ScoreType getPathScore(
    const AlignmentScores<ScoreType>& scores, const ALIGNPATH::path_t& apath, const bool isScoreOffEdge)
{
  using namespace ALIGNPATH;

  ScoreType val(0);

  // note that the intent of this function was to replicate the underling aligner and thus
  // not penalize insert-delete transitions, however as written the open score is added twice for
  // such an event. This turns out to perform much better for the performance of this tool as
  // a variant 'arm' validator. so the unintended code is staying in place for now.
  //
  /// TODO: reevaluate policy for insertion-deletion state transition score

#ifdef DEBUG_PATHSCORE
  log_os << __FUNCTION__ << " path: " << apath << "\n";
#endif

  for (const path_segment& ps : apath) {
    bool isIndel(false);  // placement of isIndel inside of this loop is the 'bug'
    switch (ps.type) {
    case MATCH:
      // if MATCH segments exist, then you're using the wrong type of CIGAR for this function
      assert(false && "Unexpected MATCH segment");
      break;
    case SEQ_MATCH:
      val += (scores.match * ps.length);
      isIndel = false;
      break;
    case SEQ_MISMATCH:
      val += (scores.mismatch * ps.length);
      isIndel = false;
      break;
    case INSERT:
    case DELETE:
      if (!isIndel) val += scores.open;
      val += (scores.extend * ps.length);
      isIndel = true;
      break;
    case SOFT_CLIP:
      if (isScoreOffEdge) val += (scores.offEdge * ps.length);
      isIndel = false;
      break;
    default:
      break;
    }

#ifdef DEBUG_PATHSCORE
    log_os << __FUNCTION__ << " path.type=" << ps.type << " path.length=" << ps.length << " val=" << val
           << "\n";
#endif
  }
  return val;
}

template <typename ScoreType>
ScoreType getMaxPathScore(
    const AlignmentScores<ScoreType>& scores,
    const ALIGNPATH::path_t&          apath,
    unsigned&                         maxReadOffset,
    unsigned&                         maxRefOffset,
    const bool                        isScoreOffEdge)
{
  using namespace ALIGNPATH;

  ScoreType val(0);
  unsigned  readOffset(0);
  unsigned  refOffset(0);

  ScoreType maxVal(0);
  maxReadOffset = 0;
  maxRefOffset  = 0;

  for (const path_segment& ps : apath) {
    bool isIndel(false);  // unintended 'bug' with positive results, see TODO note above
    switch (ps.type) {
    case MATCH:
      assert(false && "Unexpected MATCH segment");  // if MATCH segments exist, then you're using the wrong
                                                    // type of CIGAR for this function
      break;
    case SEQ_MATCH:
      val += (scores.match * ps.length);
      readOffset += ps.length;
      refOffset += ps.length;
      isIndel = false;
      break;
    case SEQ_MISMATCH:
      val += (scores.mismatch * ps.length);
      readOffset += ps.length;
      refOffset += ps.length;
      isIndel = false;
      break;
    case INSERT:
      if (!isIndel) val += scores.open;
      val += (scores.extend * ps.length);
      readOffset += ps.length;
      isIndel = true;
      break;
    case DELETE:
      if (!isIndel) val += scores.open;
      val += (scores.extend * ps.length);
      refOffset += ps.length;
      isIndel = true;
      break;
    case SOFT_CLIP:
      if (isScoreOffEdge) val += (scores.offEdge * ps.length);
      readOffset += ps.length;
      isIndel = false;
      break;
    default:
      break;
    }

    if (val > maxVal) {
      maxVal        = val;
      maxReadOffset = readOffset;
      maxRefOffset  = refOffset;
    }
  }
  return maxVal;
}
