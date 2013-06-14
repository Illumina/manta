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

#include "EstimateSVLoci.hh"
#include "ESLOptions.hh"


static
void
runESL(const ESLOptions& opt) {

    // placeholder
}



void
EstimateSVLoci::
runInternal(int argc, char* argv[]) const {

    ESLOptions opt;

    parseESLOptions(*this,argc,argv,opt);
    runESL(opt);
}
