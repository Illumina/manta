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

#include "blt_util/RegionTracker.hpp"

bool RegionTracker::isIntersectRegionImpl(const pos_t beginPos, const pos_t endPos) const
{
  if (_regions.empty()) return false;

  // 1. find first region where region.endPos > query.beginPos
  const auto posIter(_regions.upper_bound(known_pos_range2(beginPos, beginPos)));
  if (posIter == _regions.end()) return false;

  // 2. conclusion based on non-overlapping region constraint
  return (posIter->begin_pos() < endPos);
}

bool RegionTracker::isSubsetOfRegionImpl(const pos_t beginPos, const pos_t endPos) const
{
  if (_regions.empty()) return false;

  // 1. find first region where region.endPos > query.beginPos
  const auto posIter(_regions.upper_bound(known_pos_range2(beginPos, beginPos)));
  if (posIter == _regions.end()) return false;
  if (posIter->end_pos() < endPos) return false;

  // 2. conclusion based on non-overlapping region constraint
  return (posIter->begin_pos() <= beginPos);
}

void RegionTracker::addRegion(known_pos_range2 range)
{
  // check for potential set of intersecting ranges,
  // if found expand range size to represent intersection
  // remove previous content:
  const auto startOlap(_regions.upper_bound(known_pos_range2(range.begin_pos() - 1, range.begin_pos() - 1)));
  if (startOlap != _regions.end() && startOlap->begin_pos() <= (range.begin_pos() - 1)) {
    // start intersects range:
    range.set_begin_pos(startOlap->begin_pos());
  }
  auto endOlap(_regions.upper_bound(known_pos_range2(range.end_pos(), range.end_pos())));
  if (endOlap != _regions.end() && endOlap->begin_pos() <= (range.end_pos())) {
    // end intersects range:
    range.set_end_pos(endOlap->end_pos());
    endOlap++;
  }
  _regions.erase(startOlap, endOlap);
  _regions.insert(range);
}

void RegionTracker::removeToPos(const pos_t pos)
{
  auto       iter(_regions.begin());
  const auto endIter(_regions.end());
  while ((iter != endIter) && (iter->end_pos() <= (pos + 1))) {
    ++iter;
  }
  _regions.erase(_regions.begin(), iter);
}

void RegionTracker::dump(std::ostream& os) const
{
  os << "RegionTracker\n";
  for (const auto& val : _regions) {
    os << "region: " << val << "\n";
  }
}

void RegionTracker::merge(const RegionTracker& rhs)
{
  for (const auto& val : rhs._regions) {
    addRegion(val);
  }
}
