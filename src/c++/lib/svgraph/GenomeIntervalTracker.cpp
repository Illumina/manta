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

#include "GenomeIntervalTracker.hpp"

#include <iostream>

void GenomeIntervalTracker::dump(std::ostream& os) const
{
  const unsigned contigIndexSize(_regions.size());

  for (unsigned contigIndex(0); contigIndex < contigIndexSize; ++contigIndex) {
    os << "--- Regions tracked for tid " << contigIndex << "\n";
    _regions[contigIndex].dump(os);
  }
}

void GenomeIntervalTracker::merge(const GenomeIntervalTracker& rhs)
{
  const unsigned rhsContigIndexSize(rhs._regions.size());
  if (rhsContigIndexSize > _regions.size()) _regions.resize(rhsContigIndexSize);

  for (unsigned contigIndex(0); contigIndex < rhsContigIndexSize; ++contigIndex) {
    _regions[contigIndex].merge(rhs._regions[contigIndex]);
  }
}

bool GenomeIntervalTracker::isIntersectRegion(const GenomeInterval& gi) const
{
  assert(gi.tid >= 0);
  if (static_cast<unsigned>(gi.tid) >= _regions.size()) return false;
  return _regions[gi.tid].isIntersectRegion(gi.range);
}

bool GenomeIntervalTracker::isSubsetOfRegion(const GenomeInterval& gi) const
{
  assert(gi.tid >= 0);
  if (static_cast<unsigned>(gi.tid) >= _regions.size()) return false;
  return _regions[gi.tid].isSubsetOfRegion(gi.range);
}
