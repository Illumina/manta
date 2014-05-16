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

#include "GetAlignmentStats.hh"

#include "AlignmentStatsOptions.hh"

#include "blt_util/log.hh"
#include "common/OutStream.hh"
#include "manta/ReadGroupStatsUtil.hh"

#include "boost/foreach.hpp"

#include <cstdlib>



static
void
runAlignmentStats(const AlignmentStatsOptions& opt)
{
    // calculate fragment size statistics for all read groups in all bams

    // instantiate early to test for filename/permissions problems
    if (opt.alignFileOpt.alignmentFilename.empty())
    {
        log_os << "ERROR: No input files specified.\n";
        exit(EXIT_FAILURE);
    }

    const std::vector<bool>& isAlignmentTumor(opt.alignFileOpt.isAlignmentTumor);

    ReadGroupStatsSet rstats;
    const unsigned bamCount(opt.alignFileOpt.alignmentFilename.size());
    for (unsigned bamIndex(0); bamIndex<bamCount; ++bamIndex)
    {
        const std::string& bamFile(opt.alignFileOpt.alignmentFilename[bamIndex]);
        const bool isTumor(isAlignmentTumor[bamIndex]);
        extractReadGroupStatsFromBam(bamFile, rstats, isTumor);
    }

    rstats.save(opt.outputFilename.c_str());
}


void
GetAlignmentStats::
runInternal(int argc, char* argv[]) const
{
    AlignmentStatsOptions opt;

    parseAlignmentStatsOptions(*this,argc,argv,opt);
    runAlignmentStats(opt);
}
