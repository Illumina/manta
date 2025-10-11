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
/// \brief Align a contig across two breakend regions
///

#pragma once

#include "Alignment.hpp"
#include "blt_util/align_path.hpp"

struct AlignerUtil {
  static void updatePath(ALIGNPATH::path_t& path, ALIGNPATH::path_segment& ps, ALIGNPATH::align_t atype)
  {
    if (ps.type == atype) return;
    if (ps.type != ALIGNPATH::NONE) path.push_back(ps);
    ps.type   = atype;
    ps.length = 0;
  }
};

/// bookkeeping variables used during alignment backtrace
template <typename ScoreType>
struct BackTrace {
  ScoreType           max        = 0;
  AlignState::index_t state      = AlignState::MATCH;
  unsigned            queryBegin = 0;
  unsigned            refBegin   = 0;
  bool                isInit     = false;
};

/// track values needed to run the alignment backtrace:
template <typename ScoreType>
void updateBacktrace(
    const ScoreType           thisMax,
    const unsigned            refIndex,
    const unsigned            queryIndex,
    BackTrace<ScoreType>&     btrace,
    const AlignState::index_t state = AlignState::MATCH)
{
  if ((!btrace.isInit) || (thisMax > btrace.max)) {
    btrace.max        = thisMax;
    btrace.refBegin   = refIndex;
    btrace.queryBegin = queryIndex;
    btrace.isInit     = true;
    btrace.state      = state;
  }
}
