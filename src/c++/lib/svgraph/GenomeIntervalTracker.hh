// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Manta
// Copyright (c) 2013-2014 Illumina, Inc.
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

#include "GenomeInterval.hh"
#include "blt_util/RegionTracker.hh"

#include <iosfwd>


struct GenomeIntervalTracker
{
    void
    clear()
    {
        for (auto& r : _regions)
        {
            r.clear();
        }
    }

    void
    addInterval(
        const GenomeInterval& gi)
    {
        assert(gi.tid >= 0);
        if (static_cast<unsigned>(gi.tid) >= _regions.size()) _regions.resize(gi.tid+1);
        _regions[gi.tid].addRegion(gi.range);
    }

    bool
    isSubsetOfRegion(
        const GenomeInterval& gi) const
    {
        assert(gi.tid >= 0);
        if (static_cast<unsigned>(gi.tid) >= _regions.size()) return false;
        return _regions[gi.tid].isSubsetOfRegion(gi.range);
    }

private:
    std::vector<RegionTracker> _regions;
};
