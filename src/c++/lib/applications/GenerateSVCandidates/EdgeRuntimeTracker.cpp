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

#include "EdgeRuntimeTracker.hpp"

#include <iomanip>
#include <iostream>
#include <sstream>

EdgeRuntimeTracker::EdgeRuntimeTracker(std::shared_ptr<SynchronizedOutputStream> streamPtr)
  : _streamPtr(streamPtr),
    _candidateCount(0),
    _complexCandidateCount(0),
    _assembledCandidateCount(0),
    _assembledComplexCandidateCount(0)
{
}

void EdgeRuntimeTracker::stop(const EdgeInfo& edge)
{
  _edgeTime.stop();
  if (!_streamPtr) return;

  const double lastTime(_edgeTime.getWallSeconds());

  /// the purpose of the log is to identify the most troublesome cases only, so cutoff the output at a minimum
  /// time:
  static const double minLogTime(0.5);
  if (lastTime >= minLogTime) {
    std::ostringstream oss;
    oss << std::setprecision(4);
    edge.write(oss);
    oss << '\t' << lastTime;
    oss << '\t' << _candidateCount;
    oss << '\t' << _complexCandidateCount;
    oss << '\t' << _assembledCandidateCount;
    oss << '\t' << _assembledComplexCandidateCount;
    oss << '\t' << candidacyTime.getWallSeconds();
    oss << '\t' << assemblyTime.getWallSeconds();
    oss << '\t' << remoteReadRetrievalTime.getWallSeconds();
    oss << '\t' << scoreTime.getWallSeconds();
    oss << '\n';
    _streamPtr->write(oss.str());
  }
}
