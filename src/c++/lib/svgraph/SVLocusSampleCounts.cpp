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

#include "blt_util/io_util.hh"
#include "svgraph/SVLocusSampleCounts.hh"

#include <iomanip>
#include <iostream>



static
void
writeLine(
    std::ostream& os,
    const char* label,
    const double val,
    const double total)
{
    static const char sep('\t');

    os << std::fixed;
    os << label << sep;
    os << std::setprecision(0);
    os << val << sep;
    os << std::setprecision(4);
    os << val/total << '\n';
}



void
SampleReadInputCounts::
write(
    std::ostream& os) const
{
    const double dtotal(total());
    StreamScoper ss(os);
    writeLine(os,"MinMapq",minMapq,dtotal);
    writeLine(os,"Anomalous",anom,dtotal);
    writeLine(os,"AssemblyEvidence",assm,dtotal);
    writeLine(os,"Ignored",nonAnom,dtotal);
    writeLine(os,"AnomalousRemotes",remoteRecoveryCandidates,dtotal);
}



void
SampleEvidenceCounts::
write(
    std::ostream& os) const
{
    static const char sep('\t');

    double total(0);
    for (unsigned i(0); i<SVEvidenceType::SIZE; ++i)
    {
        total += eType[i];
    }

    StreamScoper ss(os);
    os << std::fixed << std::setprecision(4);
    for (unsigned i(0); i<SVEvidenceType::SIZE; ++i)
    {
        os << "EvidenceType_" << SVEvidenceType::label(i) << sep << eType[i] << sep << eType[i]/total << '\n';
    }
    os << "ClosePairs" << sep << closeCount << '\n';
}


void
SampleCounts::
write(
    std::ostream& os,
    const char* label) const
{
    os << "\n[" << label << "]\n";
    input.write(os);
    evidence.write(os);
}

