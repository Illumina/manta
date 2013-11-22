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

#include "svgraph/EdgeInfo.hh"

#include "boost/utility.hpp"

#include <ctime>

#include <iosfwd>


/// simple edge time tracker and reporter
struct EdgeRuntimeTracker : private boost::noncopyable
{
    EdgeRuntimeTracker(const std::string& outputFile);

    ~EdgeRuntimeTracker();

    void
    start(const EdgeInfo& edge);

    void
    stop();

    double
    getLastEdgeTime() const
    {
        return _lastTime;
    }

private:
    std::ostream* _osPtr;
    bool _isStart;
    clock_t _startTime;
    double _lastTime;
};
