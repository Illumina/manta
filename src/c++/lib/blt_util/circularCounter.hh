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
// <https://github.com/sequencing/licenses/>
//

///
/// \author Chris Saunders
///

#pragma once

#include <cassert>

#include <vector>


/// a fixed size circular buffer which can have true or
/// false pushed into it and the total true count
/// queried at any point
///
struct circularCounter
{
    circularCounter(
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
            assert(_count>0);
            _count--;
        }
        if (val) _count++;
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
