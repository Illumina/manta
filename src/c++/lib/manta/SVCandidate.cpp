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

#include "manta/SVCandidate.hpp"

#include <iostream>

std::ostream& operator<<(std::ostream& os, const SVCandidate& svc)
{
  static const char indent('\t');
  os << "SVCandidate:\n"
     << indent << "isImprecise?: " << svc.isImprecise() << "\n"
     << indent << "forwardTranscriptStrandReadCount: " << svc.forwardTranscriptStrandReadCount
     << " ; reverseTranscriptStrandReadCount: " << svc.reverseTranscriptStrandReadCount << "\n"
     << indent << "index candidate:assemblyAlign:assemblySegment: " << svc.candidateIndex << ":"
     << svc.assemblyAlignIndex << ":" << svc.assemblySegmentIndex << "\n";
  if (!svc.isImprecise()) {
    os << indent << "Alignment: " << svc.insertAlignment << "\n"
       << indent << "BreakendInsertSeq: " << svc.insertSeq << "\n";
  }
  if (svc.isUnknownSizeInsertion) {
    os << indent << "UnknownSizeInsertLeftSide: " << svc.unknownSizeInsertionLeftSeq << "\n"
       << indent << "UnknownSizeInsertRightSide: " << svc.unknownSizeInsertionRightSeq << "\n";
  }
  os << indent << svc.bp1 << "\n" << indent << svc.bp2 << "\n";

  return os;
}

std::ostream& operator<<(std::ostream& os, const SVObservation& svc)
{
  os << static_cast<SVCandidate>(svc);
  os << "SVObservation etype: " << SVEvidenceType::label(svc.svEvidenceType)
     << " fragtype: " << SourceOfSVEvidenceInDNAFragment::label(svc.dnaFragmentSVEvidenceSource) << "\n";
  return os;
}
