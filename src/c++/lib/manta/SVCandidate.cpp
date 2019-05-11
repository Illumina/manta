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
