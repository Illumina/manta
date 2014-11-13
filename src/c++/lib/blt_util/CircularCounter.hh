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

#include <cassert>

#include <vector>


/// A circular buffer of fixed size, S
///
/// - true/false values can be pushed in
/// - total true count among the last S pushes can be queried at any point
///    - count() is O(1) operation
///
struct CircularCounter
{
    CircularCounter(
        const unsigned initSize) :
        _count(0),
        _headPos(0),
        _data(initSize,false)
    {}

    void
    push(const bool val)
    {
        if (_data[_headPos])
        {
            if (!val)
            {
                assert(_count>0);
                _count--;
            }
        }
        else
        {
            if (val) _count++;
        }
        _data[_headPos] = val;
        _headPos = nextPos();
    }

    unsigned
    count() const
    {
        return _count;
    }

    unsigned
    size() const
    {
        return _data.size();
    }

private:

    unsigned
    nextPos() const
    {
        const unsigned pos(_headPos+1);
        if (pos>=size()) return 0;
        return pos;
    }

    unsigned _count;
    unsigned _headPos;
    std::vector<bool> _data;
};
