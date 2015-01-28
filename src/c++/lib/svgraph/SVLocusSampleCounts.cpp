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
    const char* label1,
    const char* label2,
    const double val,
    const double total)
{
    static const char sep('\t');

    os << std::fixed;
    os << label1 << '_' << label2 << ':' << sep;
    os << std::setprecision(0);
    os << val << sep;
    os << std::setprecision(4);
    os << val/total << '\n';
}



void
SampleReadInputCounts::
write(
    std::ostream& os,
    const char* label) const
{
    const double dtotal(total());
    StreamScoper ss(os);
    writeLine(os,label,"minMapq",minMapq,dtotal);
    writeLine(os,label,"Anomalous",anom,dtotal);
    writeLine(os,label,"AssemblyEvidence",assm,dtotal);
    writeLine(os,label,"Ignored",nonAnom,dtotal);
    writeLine(os,label,"AnomalousRemotes",remoteRecoveryCandidates,dtotal);
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

    StreamScoper ss(os);
    os << std::fixed << std::setprecision(4);
    for (unsigned i(0); i<SVEvidenceType::SIZE; ++i)
    {
        os << label << "_EvidenceType_" << SVEvidenceType::label(i) << ':' << sep << eType[i] << sep << eType[i]/total << '\n';
    }
    os << label << "_closePairs:" << sep << closeCount << '\n';
}
