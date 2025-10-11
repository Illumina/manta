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
/// \author Chris Saunders and Xiaoyu Chen
///

#pragma once

#include "blt_util/blt_types.hpp"

/// shared options related to read pair support:
struct PairOptions {
  explicit PairOptions(const bool isRNA)
    : minFragSupport(50), minFragProb(0.0001f), RNA(isRNA), useProperPairFlag(isRNA)

  {
  }

  /// we're interested in any fragments which cross center pos with at least N bases of support on each side
  /// (note this definition is certain to overlap the split read definition whenever N is less than the read
  /// length
  ///
  /// for reads shorter than this length, the whole read is required...
  const pos_t minFragSupport;

  const float minFragProb;

  /// This is an RNA dataset. Insert size distribution is less meaningful due to splicing
  const bool RNA;

  /// Use the proper pair flag to decide if read pairs are are anomalous for ref and alt scoring
  /// Used for RNA, since splicing makes the pair insert size distribution less informative
  const bool useProperPairFlag;
};
