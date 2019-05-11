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
/// \author Naoki Nariai
///

#include "manta/SVCandidateUtil.hpp"

bool isSVBelowMinSize(const SVCandidate& sv, const unsigned minSize)
{
  if (sv.bp1.interval.tid != sv.bp2.interval.tid) return false;

  // check if the rare IMPRECISE case arising when CIEND is a subset of CIPOS,
  // which is considered as small SV size
  if (isInvalidBreakpointInterval(sv)) return true;

  const pos_t bpSize(std::abs(sv.bp1.interval.range.center_pos() - sv.bp2.interval.range.center_pos()) - 1);
  const pos_t insertSize(sv.insertSeq.size());

  return (std::max(bpSize, insertSize) < static_cast<pos_t>(minSize));
}

bool isInvalidBreakpointInterval(const SVCandidate& sv)
{
  using namespace EXTENDED_SV_TYPE;
  const index_t svType(getExtendedSVType(sv));

  // if SV is not translocation, and IMPRECISE, then check if the SV's CIEND is a subset of CIPOS
  if (!isSVTransloc(svType) && sv.isImprecise()) {
    const bool isBp1First(sv.bp1.interval.range.begin_pos() <= sv.bp2.interval.range.begin_pos());

    const SVBreakend& bpA(isBp1First ? sv.bp1 : sv.bp2);
    const SVBreakend& bpB(isBp1First ? sv.bp2 : sv.bp1);

    return bpB.interval.range.center_pos() <= bpA.interval.range.center_pos();
  }

  return false;
}

bool isCis(const SVCandidate& sv)
{
  if (sv.bp1.interval.tid != sv.bp2.interval.tid) return false;
  if (isSameOrientation(sv.bp1.state, sv.bp2.state)) return false;
  const bool bp1Left(sv.bp1.interval.range.center_pos() < sv.bp2.interval.range.center_pos());
  if ((sv.bp1.state == SVBreakendState::RIGHT_OPEN) && bp1Left) return true;
  if ((sv.bp1.state == SVBreakendState::LEFT_OPEN) && !bp1Left) return true;
  return false;
}

SV_TYPE::index_t getSVType(const SVCandidate& sv)
{
  using namespace SV_TYPE;

  // remove failed local assemblies first:
  if ((sv.bp1.state == SVBreakendState::UNKNOWN) || (sv.bp2.state == SVBreakendState::UNKNOWN)) {
    return UNKNOWN;
  }

  const bool isBp1First(sv.bp1.interval.range.begin_pos() <= sv.bp2.interval.range.begin_pos());
  const bool isBp2First(sv.bp2.interval.range.begin_pos() <= sv.bp1.interval.range.begin_pos());

  if (sv.bp1.interval.tid != sv.bp2.interval.tid) {
    return INTERTRANSLOC;
  } else if (SVBreakendState::isSameOrientation(sv.bp1.state, sv.bp2.state)) {
    return INVERSION;
  } else if (isBp1First || isBp2First) {
    if (isInnies(isBp1First, sv.bp1.state, sv.bp2.state)) {
      return INDEL;
    } else if (isOutties(isBp1First, sv.bp1.state, sv.bp2.state)) {
      return TANDUP;
    }
  }

  return UNKNOWN;
}

namespace EXTENDED_SV_TYPE {

/// is an indel classified as insert or delete?
static index_t classifyIndel(const SVCandidate& sv)
{
  if (sv.isUnknownSizeInsertion) return INSERT;

  const bool isBp1First(sv.bp1.interval.range.begin_pos() <= sv.bp2.interval.range.begin_pos());

  const SVBreakend& bpA(isBp1First ? sv.bp1 : sv.bp2);
  const SVBreakend& bpB(isBp1First ? sv.bp2 : sv.bp1);

  const unsigned deleteSize(bpB.interval.range.begin_pos() - bpA.interval.range.begin_pos());
  const unsigned insertSize(sv.insertSeq.size());

  return ((deleteSize >= insertSize) ? DELETE : INSERT);
}

}  // namespace EXTENDED_SV_TYPE

EXTENDED_SV_TYPE::index_t getExtendedSVType(const SVCandidate& sv, const bool isForceIntraChromBnd)
{
  using namespace EXTENDED_SV_TYPE;

  const SV_TYPE::index_t svType(getSVType(sv));

  if (svType == SV_TYPE::INTERTRANSLOC) return INTERTRANSLOC;

  if (isForceIntraChromBnd) return INTRATRANSLOC;

  switch (svType) {
  case SV_TYPE::INVERSION:
    return INVERSION;
  case SV_TYPE::TANDUP:
    return TANDUP;
  case SV_TYPE::INDEL:
    return classifyIndel(sv);
  default:
    return UNKNOWN;
  }
}
