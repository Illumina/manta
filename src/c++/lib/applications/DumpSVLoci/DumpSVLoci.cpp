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

#include "manta/SVLocusSet.hh"

#include <iostream>



static
void
runDSL(const DSLOptions& opt) {

    SVLocusSet set;

    set.load(opt.graphFilename.c_str());
    set.dump(std::cout);
}



void
DumpSVLoci::
runInternal(int argc, char* argv[]) const {

    DSLOptions opt;

    parseDSLOptions(*this,argc,argv,opt);
    runDSL(opt);
}
