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
    writeLine(os,"MinMapqFiltered",minMapq,dtotal);
    writeLine(os,"NotFiltered",evidenceCount.total,dtotal);
    writeLine(os,"NotFilteredAndIgnored",evidenceCount.ignored,dtotal);
    writeLine(os,"NotFilteredAndAnomalousPair",evidenceCount.anom,dtotal);
    writeLine(os,"NotFilteredAndAnomalousPairRemotes",evidenceCount.remoteRecoveryCandidates,dtotal);
    writeLine(os,"NotFilteredAndSplitRead",evidenceCount.split,dtotal);
    writeLine(os,"NotFilteredAndLargeIndel",evidenceCount.indel,dtotal);
    writeLine(os,"NotFilteredAndSemiAligned",evidenceCount.assm,dtotal);
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

