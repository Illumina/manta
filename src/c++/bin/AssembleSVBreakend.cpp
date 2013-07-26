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

#include "applications/AssembleSVBreakend/AssembleSVBreakend.hh"


int
main(int argc, char* argv[])
{
    return AssembleSVBreakend().run(argc,argv);
    return 0;
}
