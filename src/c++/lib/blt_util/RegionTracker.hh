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

#pragma once

#include "blt_util/known_pos_range2.hh"

#include "boost/optional.hpp"

#include <iosfwd>
#include <map>
#include <set>


/// sort pos range using end_pos as the primary sort key
struct PosRangeEndSort
{
    bool
    operator()(
        const known_pos_range2& lhs,
        const known_pos_range2& rhs) const
    {
        if (lhs.end_pos() < rhs.end_pos()) return true;
        if (lhs.end_pos() == rhs.end_pos())
        {
            if (lhs.begin_pos() < rhs.begin_pos()) return true;
        }
        return false;
    }
};


/// facilitate 'rolling' region tracking and position intersect queries
///
struct RegionTracker
{
    bool
    empty() const
    {
        return _regions.empty();
    }

    void
    clear()
    {
        _regions.clear();
    }

    /// is single position in a tracked region?
    bool
    isIntersectRegion(
        const unsigned pos) const
    {
        return isIntersectRegionImpl(pos,pos+1);
    }

    /// does range intersect any tracked region?
    bool
    isIntersectRegion(
        const known_pos_range2 range) const
    {
        return isIntersectRegionImpl(range.begin_pos(),range.end_pos());
    }

    /// is range entirely contained in a region?
    bool
    isSubsetOfRegion(
        const known_pos_range2 range) const
    {
        return isSubsetOfRegionImpl(range.begin_pos(),range.end_pos());
    }

    /// add region
    ///
    /// any overlaps and adjacencies with existing regions in the tracker will be collapsed
    void
    addRegion(
        known_pos_range2 range);

    /// remove all regions which end (inclusive) before pos+1
    void
    removeToPos(
        const unsigned pos);

    // debug util
    void
    dump(std::ostream& os) const;

    typedef std::set<known_pos_range2,PosRangeEndSort>  region_t;

private:

    bool
    isIntersectRegionImpl(
        const pos_t beginPos,
        const pos_t endPos) const;

    bool
    isSubsetOfRegionImpl(
        const pos_t beginPos,
        const pos_t endPos) const;

    region_t _regions;
};


/// facilitate 'rolling' region tracking and position intersect queries
///
/// this version of RegionTracker carries a payload associated with each region
///
template <typename T>
struct RegionPayloadTracker
{
    bool
    empty() const
    {
        return _regions.empty();
    }

    boost::optional<T>
    isPayloadInRegion(const unsigned pos) const;

    /// add region
    ///
    /// any non-conflicting overlaps and adjacencies with existing regions in the tracker will be collapsed
    ///
    /// \returns false when there is an overlapping payload conflict. in this case the region is not inserted
    bool
    addRegion(
        known_pos_range2 range,
        const T payload);

    /// remove all regions which end (inclusive) before pos+1
    void
    removeToPos(
        const unsigned pos);

    // debug util
    void
    dump(std::ostream& os) const;

    typedef typename std::map<known_pos_range2,T,PosRangeEndSort> region_t;

private:
    region_t _regions;
};


#include "RegionTrackerImpl.hh"

