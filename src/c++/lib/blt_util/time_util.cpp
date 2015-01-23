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

#include "time_util.hh"

#include <cassert>



void
TimeTracker::
start()
{
    //assert((! _isStart) && "clock is running");
    _isStart = true;
    _startTime = clock();
}



void
TimeTracker::
stop()
{
    //assert(_isStart && "no clock to stop");

    static const double clockFactor(1./static_cast<double>(CLOCKS_PER_SEC));
    _totalSecs += (( clock() - _startTime )*clockFactor);
    _isStart = false;
}



double
TimeTracker::
getSeconds() const
{
    //assert((! _isStart) && "clock is running");
    return _totalSecs;
}
