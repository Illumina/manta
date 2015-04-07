// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Manta
// Copyright (c) 2013-2015 Illumina, Inc.
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

#include "FatSVCandidate.hh"

#include <iostream>



std::ostream&
operator<<(
    std::ostream& os,
    const FatSVCandidate& svc)
{
    os << static_cast<SVCandidate>(svc);
    for (unsigned eIndex(0); eIndex<SVEvidenceType::SIZE; ++eIndex)
    {
        os << "Index count for Etype: " << SVEvidenceType::label(eIndex)
           << " bp1: " << svc.bp1EvidenceIndex[eIndex].size()
           << " bp2: " << svc.bp2EvidenceIndex[eIndex].size()
           << "\n";
    }
    return os;
}
