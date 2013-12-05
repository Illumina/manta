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

#include "EdgeRuntimeTracker.hh"

#include "common/Exceptions.hh"

#include "boost/foreach.hpp"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>



EdgeRuntimeTracker::
EdgeRuntimeTracker(const std::string& outputFile) :
    _osPtr(NULL),
    _isStart(false),
    _startTime(0),
    _lastTime(0.),
    _cand(0),
    _compCand(0),
    _assmCand(0),
    _assmCompCand(0)
{
    if (outputFile.empty()) return;
    _osPtr = new std::ofstream(outputFile.c_str());
    if (! *_osPtr)
    {
        std::ostringstream oss;
        oss << "ERROR: Can't open output file: " << outputFile << '\n';
        BOOST_THROW_EXCEPTION(illumina::common::LogicException(oss.str()));
    }

    *_osPtr << std::setprecision(4);
}



EdgeRuntimeTracker::
~EdgeRuntimeTracker()
{
    if (NULL != _osPtr) delete _osPtr;
}



void
EdgeRuntimeTracker::
stop(const EdgeInfo& edge)
{
    if (! _isStart) return;

    static const double clockFactor(1./static_cast<double>(CLOCKS_PER_SEC));
    _lastTime = (( clock() - _startTime )*clockFactor);
    _isStart = false;

    /// the purpose of the log is to identify the most troublesome cases only, so cutoff the output at a minimum time:
    static const double minLogTime(0.5);
    if (_lastTime >= minLogTime)
    {
        if (NULL != _osPtr)
        {
            edge.write(*_osPtr);
            *_osPtr << '\t' << _lastTime
                    << '\t' << _cand
                    << '\t' << _compCand
                    << '\t' << _assmCand
                    << '\t' << _assmCompCand
                    << '\n';
        }
    }
}
