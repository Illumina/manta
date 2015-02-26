// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Starka
// Copyright (c) 2009-2014 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/sequencing/licenses/>
//

///
/// \author Chris Saunders
///

#include "blt_util/RegionTracker.hh"



bool
RegionTracker::
isIntersectRegionImpl(
    const pos_t beginPos,
    const pos_t endPos) const
{
    if (_regions.empty()) return false;

    // 1. find first region where region.endPos > query.beginPos
    const auto posIter(_regions.upper_bound(known_pos_range2(beginPos,beginPos)));
    if (posIter == _regions.end()) return false;

    // 2. conclusion based on non-overlapping region constraint
    return (posIter->begin_pos() < endPos);
}


bool
RegionTracker::
isSubsetOfRegionImpl(
    const pos_t beginPos,
    const pos_t endPos) const
{
    if (_regions.empty()) return false;

    // 1. find first region where region.endPos > query.beginPos
    const auto posIter(_regions.upper_bound(known_pos_range2(beginPos,beginPos)));
    if (posIter == _regions.end()) return false;
    if (posIter->end_pos() < endPos) return false;

    // 2. conclusion based on non-overlapping region constraint
    return (posIter->begin_pos() <= beginPos);
}



void
RegionTracker::
addRegion(known_pos_range2 range)
{
    // check for potential set of intersecting ranges,
    // if found expand range size to represent intersection
    // remove previous content:
    const auto startOlap(_regions.upper_bound(known_pos_range2(range.begin_pos()-1,range.begin_pos()-1)));
    if (startOlap != _regions.end() && startOlap->begin_pos() <= (range.begin_pos()-1))
    {
        // start intersects range:
        range.set_begin_pos(startOlap->begin_pos());
    }
    auto endOlap(_regions.upper_bound(known_pos_range2(range.end_pos(),range.end_pos())));
    if (endOlap != _regions.end() && endOlap->begin_pos() <= (range.end_pos()))
    {
        // end intersects range:
        range.set_end_pos(endOlap->end_pos());
        endOlap++;
    }
    _regions.erase(startOlap,endOlap);
    _regions.insert(range);
}



void
RegionTracker::
removeToPos(const unsigned pos)
{
    for (auto iter(_regions.begin()) ; iter != _regions.end() ; ++iter)
    {
        if (iter->end_pos() > (pos+1)) return;
        _regions.erase(iter);
    }
}



void
RegionTracker::
dump(std::ostream& os) const
{
    os << "RegionTracker\n";
    for (const auto& val : _regions)
    {
        os << "region: " << val << "\n";
    }
}
