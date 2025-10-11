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

#pragma once

#include "manta/SVCandidate.hpp"

/// SVComplexCandidate represents an associated grouping of multiple breakend pairs
///
/// Examples: The two breakend pairs of a simple inversion would form a complex candidates composed of the two
/// SVCandidates for the forward and reverse junctions.
///
/// for a complex candidates we want to test the concept that the full set of breakends occurred together as
/// part of the same event, and thus be able to call the full event with perhaps less evidence per each single
/// breakend than would be acceptable during regular calling.
///
struct SVMultiJunctionCandidate {
  SVMultiJunctionCandidate() {}

  std::vector<SVCandidate> junction;

  /// TODO: need to design a quick data structure to iterate through complex event breakend regions and pull
  /// out the associated candidate's breakends -- is this just a mapping from MJ breakend 'groups' to literal
  /// junction breakends??
};

std::ostream& operator<<(std::ostream& os, const SVMultiJunctionCandidate& svc);
