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

#include "blt_util/time_util.hh"
#include "svgraph/EdgeInfo.hh"

#include "boost/utility.hpp"


#include <iosfwd>



/// simple edge time tracker and reporter
struct EdgeRuntimeTracker : private boost::noncopyable
{
    EdgeRuntimeTracker(const std::string& outputFile);

    ~EdgeRuntimeTracker();

    void
    start()
    {
        edgeTime.clear();
        candTime.clear();
        assmTime.clear();
        scoreTime.clear();
        remoteTime.clear();

        edgeTime.resume();
        _cand = 0;
        _compCand = 0;
        _assmCand = 0;
        _assmCompCand = 0;
    }

    void
    stop(const EdgeInfo& edge);

    double
    getLastEdgeTime() const
    {
        return edgeTime.getSeconds();
    }

    void
    addCand(const bool isComplex)
    {
        if (isComplex) _compCand++;
        else           _cand++;
    }

    void
    addAssm(const bool isComplex)
    {
        if (isComplex) _assmCompCand++;
        else           _assmCand++;
    }

    TimeTracker candTime;
    TimeTracker assmTime;
    TimeTracker scoreTime;
    TimeTracker remoteTime;
private:
    std::ostream* _osPtr;
    TimeTracker edgeTime;

    unsigned _cand;
    unsigned _compCand;
    unsigned _assmCand;
    unsigned _assmCompCand;
};
