//
// Manta - Structural Variant and Indel Caller
// Copyright (c) 2013-2019 Illumina, Inc.
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

#include <iosfwd>

#include "boost/noncopyable.hpp"

#include "blt_util/io_util.hpp"
#include "blt_util/time_util.hpp"
#include "svgraph/EdgeInfo.hpp"

/// Simple edge time tracker and reporter
///
/// When multiple trackers are used across multiple threads, their output is coordinated at the file output
/// level by a shared SynchronizedOutputStream, which must be provided at construction.
///
struct EdgeRuntimeTracker : private boost::noncopyable {
  /// Default construction method to synchronize Edge tracking output over multiple threads
  explicit EdgeRuntimeTracker(std::shared_ptr<SynchronizedOutputStream> streamPtr);

  /// Simpler construction method for single-thread use
  explicit EdgeRuntimeTracker(const std::string& outputFilename)
    : EdgeRuntimeTracker(std::make_shared<SynchronizedOutputStream>(outputFilename))
  {
  }

  void start()
  {
    _edgeTime.clear();
    candidacyTime.clear();
    assemblyTime.clear();
    scoreTime.clear();
    remoteReadRetrievalTime.clear();

    _edgeTime.resume();
    _candidateCount                 = 0;
    _complexCandidateCount          = 0;
    _assembledCandidateCount        = 0;
    _assembledComplexCandidateCount = 0;
  }

  void stop(const EdgeInfo& edge);

  CpuTimes getLastEdgeTime() const { return _edgeTime.getTimes(); }

  void addCandidate(const bool isComplex)
  {
    if (isComplex)
      _complexCandidateCount++;
    else
      _candidateCount++;
  }

  void addAssembledCandidate(const bool isComplex)
  {
    if (isComplex)
      _assembledComplexCandidateCount++;
    else
      _assembledCandidateCount++;
  }

  TimeTracker candidacyTime;
  TimeTracker assemblyTime;
  TimeTracker scoreTime;

  /// Time to retrieve reads from remote locations prior to assembling large insertions.
  /// Note this is a subset of assemblyTime
  TimeTracker remoteReadRetrievalTime;

  /// TestEdgeRuntimeTracker is a friend structure of EdgeRuntimeTracker. So that it can access private
  /// members of EdgeRuntimeTracker.
  friend struct TestEdgeRuntimeTracker;

private:
  std::shared_ptr<SynchronizedOutputStream> _streamPtr;
  TimeTracker                               _edgeTime;

  unsigned _candidateCount;
  unsigned _complexCandidateCount;
  unsigned _assembledCandidateCount;
  unsigned _assembledComplexCandidateCount;
};
