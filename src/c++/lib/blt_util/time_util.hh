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

#include "boost/timer/timer.hpp"
#include "boost/utility.hpp"


/// helper functions for boost::timer::cpu_times
inline
void
merge(
    boost::timer::cpu_times& lhs,
    const boost::timer::cpu_times& rhs)
{
    lhs.wall += rhs.wall;
    lhs.user += lhs.user;
    lhs.system += lhs.system;
}

inline
void
difference(
    boost::timer::cpu_times& lhs,
    const boost::timer::cpu_times& rhs)
{
    lhs.wall -= rhs.wall;
    lhs.user -= lhs.user;
    lhs.system -= lhs.system;
}


namespace boost {
namespace serialization {

template<class Archive>
void serialize(Archive & ar, boost::timer::cpu_times & t, const unsigned /*version*/)
{
    ar & t.wall & t.user & t.system;
}

}
}


/// simple time track utility
struct TimeTracker
{
    TimeTracker() { _timer.stop(); }

    void
    clear()
    {
        _isReset = true;
    }

    /// starts clock without reset to accumulate total time
    void
    resume()
    {
        //assert((! _isStart) && "clock is running");
        if (_isReset)
        {
            _timer.start();
            _isReset = false;
        }
        else _timer.resume();
    }

    /// stop clock
    void
    stop()
    {
        _timer.stop();
    }

    boost::timer::cpu_times
    getTimes() const
    {
        return _timer.elapsed();
    }

    /// DEPRECATED get user cpu time in seconds
    ///
    /// timer must be stopped
    double
    getSeconds() const
    {
        using namespace boost::chrono;
        return static_cast<double>(duration_cast<milliseconds>(nanoseconds(getTimes().user)).count())/1000.;
    }

private:
    bool _isReset = true;
    boost::timer::cpu_timer _timer;
};


/// utility for timetracker for scope based start-stop scenarios:
struct TimeScoper : private boost::noncopyable
{
    TimeScoper(TimeTracker& t) : _t(t)
    {
        _t.resume();
    }

    ~TimeScoper()
    {
        _t.stop();
    }
private:
    TimeTracker& _t;
};
