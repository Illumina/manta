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

#include "manta/SVCandidate.hh"

#include <iostream>



std::ostream&
operator<<(
    std::ostream& os,
    const SVCandidate& svc)
{
    static const char indent('\t');
    os << "SVCandidate:\n"
       << indent << "isImprecise?: " << svc.isImprecise() << "\n"
       << indent << "fwReads: " << svc.fwReads << " ; rvReads: " << svc.rvReads << "\n"
       << indent << "index candidate:assemblyAlign:assemblySegment: " << svc.candidateIndex
       << ":" << svc.assemblyAlignIndex
       << ":" << svc.assemblySegmentIndex
       << "\n";
    if (! svc.isImprecise())
    {
        os << "\tAlignment: " << svc.insertAlignment << "\n"
           << "\tBreakendInsertSeq: " << svc.insertSeq << "\n";
    }
    os << "\t" << svc.bp1 << "\n"
       << "\t" << svc.bp2 << "\n";
    return os;
}



std::ostream&
operator<<(
    std::ostream& os,
    const SVObservation& svc)
{
    os << static_cast<SVCandidate>(svc);
    os << "SVObservation etype: " << SVEvidenceType::label(svc.evtype)
       << " fragtype: " << FRAGSOURCE::label(svc.fragSource) << "\n";
    return os;
}
