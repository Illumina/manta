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

#include "blt_util/blt_types.hh"
#include "blt_util/RangeMap.hh"

#include <cassert>


/// simple map of position to depth
///
/// assumes that a narrow list of positions is maintained so that
/// array based lookup optimizations can be used
///
struct depth_buffer
{
    unsigned
    val(const pos_t pos) const
    {
        return _data.getConstRefDefault(pos,0);
    }

    void
    inc(const pos_t pos)
    {
        _data.getRef(pos) += 1;
    }

    void
    clear_pos(const pos_t pos)
    {
        if (_data.isKeyPresent(pos)) _data.erase(pos);
    }

    /// return true if buffered depth exceeds depth in [begin,end]
    bool
    is_range_ge_than(const pos_t begin,
                     const pos_t end,
                     const unsigned depth) const
    {
        assert(begin <= end);
        for (pos_t i(begin); i<=end; ++i)
        {
            if (val(i) >= depth) return true;
        }
        return false;
    }

private:
    typedef RangeMap<pos_t,unsigned> count_t;
    count_t _data;
};
