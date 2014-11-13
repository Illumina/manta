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

#include "CheckSVLoci.hh"
#include "CSLOptions.hh"

#include "blt_util/log.hh"
#include "htsapi/bam_header_util.hh"
#include "svgraph/SVLocusSet.hh"

#include "boost/archive/binary_oarchive.hpp"

#include <fstream>
#include <iostream>



static
void
runCSL(const CSLOptions& opt)
{
    log_os << "INFO: loading graph: " << opt.graphFilename << "\n";

    SVLocusSet set;
    set.load(opt.graphFilename.c_str());

    log_os << "INFO: cleaning/finalizing graph: " << opt.graphFilename << "\n";

    set.finalize();

    log_os << "INFO: checking cleaned graph: " << opt.graphFilename << "\n";

    set.checkState(true,true);

    log_os << "INFO: finished checking graph: " << opt.graphFilename << "\n";
}



void
CheckSVLoci::
runInternal(int argc, char* argv[]) const
{
    CSLOptions opt;

    parseCSLOptions(*this,argc,argv,opt);
    runCSL(opt);
}
