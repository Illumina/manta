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

#include "boost/utility.hpp"

#include <ctime>



/// simple time track utility
struct TimeTracker
{
    void reset()
    {
        _isStart = false;
        _startTime = 0;
        _totalSecs = 0;
    }

    void start();

    void stop();

    /// get total time in seconds
    ///
    /// timer must be stopped
    double getSeconds() const;

private:
    bool _isStart = false;
    clock_t _startTime = 0;
    double _totalSecs = 0;
};


/// utility for timetracker for scope based start-stop scenarios:
struct TimeScoper : private boost::noncopyable
{
    TimeScoper(TimeTracker& t) : _t(t)
    {
        _t.start();
    }

    ~TimeScoper()
    {
        _t.stop();
    }
private:
    TimeTracker& _t;
};
