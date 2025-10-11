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
/// \author Naoki Nariai
///

#include "FatSVCandidate.hpp"

#include <iostream>

std::ostream& operator<<(std::ostream& os, const FatSVCandidate& svc)
{
  os << static_cast<SVCandidate>(svc);
  for (unsigned eIndex(0); eIndex < SVEvidenceType::SIZE; ++eIndex) {
    os << "Index count for Etype: " << SVEvidenceType::label(eIndex);
    for (unsigned bamIndex(0); bamIndex < svc.bp1EvidenceIndex[eIndex].size(); ++bamIndex) {
      os << "Bam index: " << bamIndex << " bp1: " << svc.bp1EvidenceIndex[eIndex][bamIndex].size()
         << " bp2: " << svc.bp2EvidenceIndex[eIndex][bamIndex].size() << "\n";
    }
  }
  return os;
}
