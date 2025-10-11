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

#include "SVEvidence.hpp"
#include "manta/SVCandidate.hpp"
#include "manta/SVScoreInfo.hpp"

/// manage all per-junction information consumed by an SV calling model
///
/// using this object facilities multi-breakend event scoring, but clearly
/// separating out per-junction input info from junction-independent info
///
struct JunctionCallInfo {
  JunctionCallInfo() : _sv(nullptr), _evidence(nullptr), _baseInfo(nullptr), _spanningPairWeight(0) {}

  const SVCandidate& getSV() const
  {
    assert(nullptr != _sv);
    return *_sv;
  }

  const SVEvidence& getEvidence() const
  {
    assert(nullptr != _evidence);
    return *_evidence;
  }

  const SVScoreInfo& getBaseInfo() const
  {
    assert(nullptr != _baseInfo);
    return *_baseInfo;
  }

  float getSpanningWeight() const { return _spanningPairWeight; }

  void init(
      const SVCandidate& sv,
      const SVEvidence&  evidence,
      const SVScoreInfo& baseInfo,
      const float        spanningPairWeight)
  {
    _sv                 = &sv;
    _evidence           = &evidence;
    _baseInfo           = &baseInfo;
    _spanningPairWeight = spanningPairWeight;
  }

private:
  const SVCandidate* _sv;
  const SVEvidence*  _evidence;
  const SVScoreInfo* _baseInfo;
  float              _spanningPairWeight;
};
