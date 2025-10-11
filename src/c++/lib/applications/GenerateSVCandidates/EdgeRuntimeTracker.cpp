//
// Manta - Structural Variant and Indel Caller
// Copyright (c) 2013-2025 Illumina, Inc.
//
// This program is licensed under the terms of the Polyform strict license
//
// ***As far as the law allows, the software comes as is, without
// any warranty or condition, and the licensor will not be liable
// to you for any damages arising out of these terms or the use
// or nature of the software, under any kind of legal claim.***
//
// You should have received a copy of the PolyForm Strict License 1.0.0
// along with this program.  If not, see <https://polyformproject.org/licenses/strict/1.0.0>.
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
