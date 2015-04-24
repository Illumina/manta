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

#include "SummarizeAlignmentStats.hh"
#include "SASOptions.hh"

#include "blt_util/log.hh"
#include "manta/ReadGroupStatsSet.hh"

#include <iostream>



static
void
runSAS(const SASOptions& opt)
{
    static const float quantLevel[] = { 0.01f, 0.05f, 0.10f, 0.25f, 0.50f, 0.75f, 0.90f, 0.95f, 0.99f };
    static const unsigned quantLevelCount(sizeof(quantLevel)/sizeof(float));

    std::ostream& report_os(std::cout);

    ReadGroupStatsSet rgss;
    rgss.load(opt.statsFilename.c_str());

    const unsigned groupCount(rgss.size());
    for (unsigned groupIndex(0); groupIndex<groupCount; ++groupIndex)
    {
        const ReadGroupStatsSet::KeyType& key(rgss.getKey(groupIndex));
#ifdef READ_GROUPS
        report_os << "bamFile:\t" << key.bamLabel << '\n';
        report_os << "readGroup:\t" << key.rgLabel << '\n';
#else
        report_os << "group:\t" << key.bamLabel << '\n';
#endif

        const ReadGroupStats& rgs(rgss.getStats(groupIndex));
        report_os << "fragment length observations:\t" << rgs.fragStats.totalObservations() << '\n';
        report_os << "fragment length quantiles:\n";
        for (unsigned quantLevelIndex(0); quantLevelIndex<quantLevelCount; ++quantLevelIndex)
        {
            report_os << quantLevel[quantLevelIndex] << '\t' << rgs.fragStats.quantile(quantLevel[quantLevelIndex]) << '\n';
        }
        report_os << '\n';
    }
}



void
SummarizeAlignmentStats::
runInternal(int argc, char* argv[]) const
{
    SASOptions opt;

    parseSASOptions(*this,argc,argv,opt);
    runSAS(opt);
}
