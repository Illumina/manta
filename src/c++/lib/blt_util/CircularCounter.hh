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
        _dataSize(0),
        _maxCount(0),
        _data(initSize,false)
    {
        assert(initSize>0);
    }

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
            if (_count > _maxCount) _maxCount = _count;
        }
        _data[_headPos] = val;
        _headPos = nextPos();
        if (_dataSize < size()) _dataSize++;
    }

    // change the value on the head of the buffer
    void
    replace(const bool val)
    {
        assert(_dataSize>0);
        _headPos=lastPos();
        _dataSize--;
        push(val);
    }

    unsigned
    count() const
    {
        return _count;
    }

    unsigned
    maxCount() const
    {
        return _maxCount;
    }

    /// less than or equal to size(), according to the number
    /// of observations pushed
    unsigned
    dataSize() const
    {
        return _dataSize;
    }

    unsigned
    size() const
    {
        return _data.size();
    }

private:
    unsigned
    lastPos() const
    {
        if (_headPos==0) return (size()-1);
        return _headPos-1;
    }

    unsigned
    nextPos() const
    {
        const unsigned pos(_headPos+1);
        if (pos>=size()) return 0;
        return pos;
    }

    unsigned _count;
    unsigned _headPos;
    unsigned _dataSize;
    unsigned _maxCount;
    std::vector<bool> _data;
};
