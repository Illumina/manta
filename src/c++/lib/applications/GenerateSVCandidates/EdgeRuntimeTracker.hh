//
// Manta - Structural Variant and Indel Caller
// Copyright (c) 2013-2018 Illumina, Inc.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//

/// \file
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
    explicit
    EdgeRuntimeTracker(const std::string& outputFile);

    ~EdgeRuntimeTracker();

    void
    start()
    {
        edgeTime.clear();
        candidacyTime.clear();
        assemblyTime.clear();
        scoreTime.clear();
        remoteReadRetrievalTime.clear();

        edgeTime.resume();
        _candidateCount = 0;
        _complexCandidateCount = 0;
        _assembledCandidateCount = 0;
        _assembledComplexCandidateCount = 0;
    }

    void
    stop(const EdgeInfo& edge);

    CpuTimes
    getLastEdgeTime() const
    {
        return edgeTime.getTimes();
    }

    void
    addCandidate(const bool isComplex)
    {
        if (isComplex) _complexCandidateCount++;
        else           _candidateCount++;
    }

    void
    addAssembledCandidate(const bool isComplex)
    {
        if (isComplex) _assembledComplexCandidateCount++;
        else           _assembledCandidateCount++;
    }

    TimeTracker candidacyTime;
    TimeTracker assemblyTime;
    TimeTracker scoreTime;

    /// Time to retrieve reads from remote locations prior to assembling large insertions.
    /// Note this is a subset of assemblyTime
    TimeTracker remoteReadRetrievalTime;

    /// TestEdgeRuntimeTracker is a friend structure of EdgeRuntimeTracker. So that it can access private
    /// members of EdgeRuntimeTracker. As we need to close the output stream which is a private member, so for
    /// unit test writing this friend structure has been created.
    friend struct TestEdgeRuntimeTracker;

private:
    std::ostream* _osPtr;
    TimeTracker edgeTime;

    unsigned _candidateCount;
    unsigned _complexCandidateCount;
    unsigned _assembledCandidateCount;
    unsigned _assembledComplexCandidateCount;
};
