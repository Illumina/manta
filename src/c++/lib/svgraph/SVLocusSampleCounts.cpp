// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Manta
// Copyright (c) 2013-2014 Illumina, Inc.
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

#include "svgraph/SVLocusSampleCounts.hh"

#include <iostream>


void
SampleReadInputCounts::
write(
    std::ostream& os,
    const char* label) const
{
    static const char sep('\t');

    const double dtotal(total());
    os << label << "_Anomalous:" << sep << anom << sep << anom/dtotal << '\n';
    os << label << "_AssemblyEvidence:" << sep << assm << sep << assm/dtotal << '\n';
    os << label << "_Ignored:" << sep << nonAnom << sep << nonAnom/dtotal << '\n';
    os << label << "_AnomalousRemotes:" << sep << remoteRecoveryCandidates << sep << remoteRecoveryCandidates/dtotal << '\n';
}



void
SampleEvidenceCounts::
write(
    std::ostream& os,
    const char* label) const
{
    static const char sep('\t');

    double total(0);
    for (unsigned i(0); i<SVEvidenceType::SIZE; ++i)
    {
        total += eType[i];
    }

    for (unsigned i(0); i<SVEvidenceType::SIZE; ++i)
    {
        os << label << "_EvidenceType_" << SVEvidenceType::label(i) << ':' << sep << eType[i] << sep << eType[i]/total << '\n';
    }
    os << label << "_closePairs:" << sep << closeCount << '\n';
}
