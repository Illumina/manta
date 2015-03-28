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

#include "CheckSVLoci.hh"
#include "CSLOptions.hh"

#include "blt_util/log.hh"
#include "svgraph/SVLocusSet.hh"



static
void
runCSL(const CSLOptions& opt)
{
    SVLocusSet set;
    set.load(opt.graphFilename.c_str());
    set.finalize();
    set.checkState(true,true);
}



void
CheckSVLoci::
runInternal(int argc, char* argv[]) const
{
    CSLOptions opt;

    parseCSLOptions(*this,argc,argv,opt);
    runCSL(opt);
}
