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

#include "EdgeRuntimeTracker.hh"

#include "common/Exceptions.hh"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>



EdgeRuntimeTracker::
EdgeRuntimeTracker(
    const std::string& outputFile) :
    _osPtr(nullptr),
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
    if (nullptr != _osPtr) delete _osPtr;
}



void
EdgeRuntimeTracker::
stop(const EdgeInfo& edge)
{
    edgeTime.stop();
    const double lastTime(edgeTime.getSeconds());

    /// the purpose of the log is to identify the most troublesome cases only, so cutoff the output at a minimum time:
    static const double minLogTime(0.5);
    if (lastTime >= minLogTime)
    {
        if (NULL != _osPtr)
        {
            edge.write(*_osPtr);
            *_osPtr << '\t' << lastTime
                    << '\t' << _cand
                    << '\t' << _compCand
                    << '\t' << _assmCand
                    << '\t' << _assmCompCand
                    << '\t' << candTime.getSeconds()
                    << '\t' << assmTime.getSeconds()
                    << '\t' << remoteTime.getSeconds()
                    << '\t' << scoreTime.getSeconds()
                    << '\n';
        }
    }
}
