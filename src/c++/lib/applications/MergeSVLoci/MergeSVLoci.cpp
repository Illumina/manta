// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Manta
// Copyright (c) 2013 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/downloads/sequencing/licenses/>.
//

///
/// \author Chris Saunders
///

#include "MergeSVLoci.hh"
#include "MSLOptions.hh"

#include "blt_util/log.hh"
#include "common/OutStream.hh"
#include "manta/SVLocusSet.hh"



static
void
runMSL(const MSLOptions& opt)
{

    {
        // early test that we have permission to write to output file
        OutStream outs(opt.outputFilename);
    }

    SVLocusSet mergedSet;

    BOOST_FOREACH(const std::string& graphFile, opt.graphFilename)
    {
        if (opt.isVerbose)
        {
            log_os << "INFO: Merging file: " << graphFile << "\n";
        }

        if (mergedSet.empty())
        {
            mergedSet.load(graphFile.c_str());
        }
        else
        {
            SVLocusSet inputSet;
            inputSet.load(graphFile.c_str());
            mergedSet.merge(inputSet);
        }

        if (opt.isVerbose)
        {
            log_os << "INFO: Finished merging file: " << graphFile << "\n";
        }
    }

    mergedSet.clean();
    mergedSet.save(opt.outputFilename.c_str());
}



void
MergeSVLoci::
runInternal(int argc, char* argv[]) const
{

    MSLOptions opt;

    parseMSLOptions(*this,argc,argv,opt);
    runMSL(opt);
}
