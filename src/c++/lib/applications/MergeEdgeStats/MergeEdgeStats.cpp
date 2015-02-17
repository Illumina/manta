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

#include "MergeEdgeStats.hh"
#include "MESOptions.hh"
#include "common/OutStream.hh"
#include "svgraph/GSCEdgeStats.hh"

#include "blt_util/log.hh"



static
void
runMES(const MESOptions& opt)
{
    {
        // early test that we have permission to write to output file(s)
        OutStream outs(opt.outputFilename);
        if (! opt.reportFilename.empty())
        {
            OutStream reps(opt.reportFilename);
        }
    }

    GSCEdgeStats mergedStats;
    bool isFirst(true);
    for (const std::string& statsFilename : opt.statsFilename)
    {
        if (isFirst)
        {
            mergedStats.load(statsFilename.c_str());
            isFirst=false;
        }
        else
        {
            GSCEdgeStats inputStats;
            inputStats.load(statsFilename.c_str());
            mergedStats.merge(inputStats);
        }

    }

    mergedStats.save(opt.outputFilename.c_str());
    if (! opt.reportFilename.empty())
    {
        mergedStats.report(opt.reportFilename.c_str());
    }
}



void
MergeEdgeStats::
runInternal(int argc, char* argv[]) const
{
    MESOptions opt;

    parseMESOptions(*this,argc,argv,opt);
    runMES(opt);
}
