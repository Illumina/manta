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

#include "SummarizeSVLoci.hh"
#include "SSLOptions.hh"

#include "svgraph/SVLocusSet.hh"

#include <iostream>



static
void
runSSL(const SSLOptions& opt)
{
    SVLocusSet set;

    set.load(opt.graphFilename.c_str());

    std::ostream& os(std::cout);

    if (opt.isGlobalStats)
    {
        set.dumpStats(os);
    }
    else
    {
        set.dumpLocusStats(os);
    }
}



void
SummarizeSVLoci::
runInternal(int argc, char* argv[]) const
{
    SSLOptions opt;

    parseSSLOptions(*this,argc,argv,opt);
    runSSL(opt);
}
