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
