// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Manta
// Copyright (c) 2013 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/downloads/sequencing/licenses/>.
//

///
/// \author Chris Saunders
///

///
/// this is the beginning of a redesign to known_pos_range
/// to be more efficient for the manta case
///

#pragma once

#include "blt_util/blt_types.hh"

#include "boost/serialization/level.hpp"

#include <algorithm>
#include <iosfwd>


/// \brief integer ranges which are right open
///
struct known_pos_range2
{
    known_pos_range2() :
        _begin_pos(0),
        _end_pos(0)
    {}

    known_pos_range2(
        const pos_t bp,
        const pos_t ep) :
        _begin_pos(bp),
        _end_pos(ep)
    {}

    void
    set_begin_pos(const pos_t pos)
    {
        _begin_pos=pos;
    }

    void
    set_end_pos(const pos_t pos)
    {
        _end_pos=pos;
    }

    void
    set_range(const pos_t begin,
              const pos_t end)
    {
        set_begin_pos(begin);
        set_end_pos(end);
    }

    pos_t
    begin_pos() const
    {
        return _begin_pos;
    }

    pos_t
    end_pos() const
    {
        return _end_pos;
    }

    pos_t
    center_pos() const
    {
        return _begin_pos + ((std::max(size(),1u)-1)/2);
    }

    bool
    is_pos_intersect(const pos_t pos) const
    {
        return ((pos >= _begin_pos) &&
                (pos < _end_pos));
    }

    bool
    is_range_intersect(const known_pos_range2& pr) const
    {
        return ((pr._end_pos > _begin_pos) &&
                (pr._begin_pos < _end_pos));
    }

    /// does this range completely overlap pr?
    bool
    is_superset_of(const known_pos_range2& pr) const
    {
        return
            ((pr._end_pos <= _end_pos) &&
             (pr._begin_pos >= _begin_pos));
    }

    unsigned
    size() const
    {
        return std::max(0,_end_pos-_begin_pos);
    }

    bool
    operator<(const known_pos_range2& rhs) const
    {
        if (_begin_pos < rhs._begin_pos) return true;
        if (_begin_pos == rhs._begin_pos)
        {
            if (_end_pos < rhs._end_pos) return true;
        }
        return false;
    }

    bool
    operator==(const known_pos_range2& rhs) const
    {
        return ((_begin_pos==rhs._begin_pos) && (_end_pos==rhs._end_pos));
    }

    // expand range to extend of a second range:
    void
    merge_range(const known_pos_range2& kpr)
    {
        if (kpr._begin_pos<_begin_pos) _begin_pos=kpr._begin_pos;
        if (kpr._end_pos>_end_pos) _end_pos=kpr._end_pos;
    }

    void
    clear()
    {
        _begin_pos=0;
        _end_pos=0;
    }


    template<class Archive>
    void serialize(Archive& ar, const unsigned /* version */)
    {
        ar& _begin_pos& _end_pos;
    }

private:
    pos_t _begin_pos;
    pos_t _end_pos;
};


std::ostream& operator<<(std::ostream& os, const known_pos_range2& pr);

BOOST_CLASS_IMPLEMENTATION(known_pos_range2, boost::serialization::object_serializable)

