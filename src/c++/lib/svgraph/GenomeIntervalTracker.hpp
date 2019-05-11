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

#pragma once

#include "GenomeInterval.hpp"
#include "blt_util/RegionTracker.hpp"

#include "boost/serialization/vector.hpp"

#include <iosfwd>

/// Accumulates genome intervals, then provides tests on how any new genome interval interacts with the
/// accumulated interval set
///
struct GenomeIntervalTracker {
  void clear() { _regions.clear(); }

  void addInterval(const GenomeInterval& gi)
  {
    assert(gi.tid >= 0);
    if (static_cast<unsigned>(gi.tid) >= _regions.size()) _regions.resize(gi.tid + 1);
    _regions[gi.tid].addRegion(gi.range);
  }

  /// Merge genome interval content from rhs into this
  void merge(const GenomeIntervalTracker& rhs);

  /// Return true if \p gi intersects any of the intervals already added to this object
  bool isIntersectRegion(const GenomeInterval& gi) const;

  /// Return true if \p gi is entirely contained within the intervals already added to this object
  bool isSubsetOfRegion(const GenomeInterval& gi) const;

  /// debug util
  void dump(std::ostream& os) const;

  template <class Archive>
  void serialize(Archive& ar, const unsigned /* version */)
  {
    ar& _regions;
  }

private:
  std::vector<RegionTracker> _regions;
};

BOOST_CLASS_IMPLEMENTATION(GenomeIntervalTracker, boost::serialization::object_serializable)
