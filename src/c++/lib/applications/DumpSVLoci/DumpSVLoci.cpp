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

#include "DumpSVLoci.hh"
#include "DSLOptions.hh"

#include "blt_util/bam_header_util.hh"
#include "svgraph/SVLocusSet.hh"

#include "boost/archive/binary_oarchive.hpp"

#include <fstream>
#include <iostream>



static
void
runDSL(const DSLOptions& opt)
{

    SVLocusSet set;
    set.load(opt.graphFilename.c_str());

    const SVLocusSet& cset(set);

    std::ostream& os(std::cout);

    // add this handy map of chromosome id to chromosome label at the start of all output types:
    os << cset.header << "\n";

    if (! opt.region.empty())
    {
        int32_t tid,beginPos,endPos;
        parse_bam_region(set.header, opt.region, tid, beginPos, endPos); // parse the region

        set.dumpRegion(os,GenomeInterval(tid,beginPos,endPos));
    }
    else if (opt.isLocusIndex)
    {
        const SVLocus& locus(cset.getLocus(opt.locusIndex));
        if (opt.locusFilename.empty())
        {
            os << locus;
        }
        else
        {
            std::ofstream ofs(opt.locusFilename.c_str(), std::ios::binary);
            boost::archive::binary_oarchive oa(ofs);
            oa << locus;
        }
    }
    else
    {
        cset.dump(os);
    }
}



void
DumpSVLoci::
runInternal(int argc, char* argv[]) const
{

    DSLOptions opt;

    parseDSLOptions(*this,argc,argv,opt);
    runDSL(opt);
}
