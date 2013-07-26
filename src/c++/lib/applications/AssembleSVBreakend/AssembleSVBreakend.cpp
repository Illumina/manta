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
/// \author Ole Schulz-Trieglaff
///

#include "AssembleSVBreakend.hh"
#include "applications/GenerateSVCandidates/GSCOptions.hh"
#include "applications/AssembleSVBreakend/ASBOptions.hh"

#include "manta/SVLocusAssembler.hh"

#include <fstream>
#include <iostream>


static
void
runASB(const ASBOptions& opt)
{
	///
    SVLocusAssembler svla(opt);
}


void
AssembleSVBreakend::
runInternal(int argc, char* argv[]) const
{
	ASBOptions opt;
    parseASBOptions(*this,argc,argv,opt);
    //parseGSCOptions(*this,argc,argv,opt);
    runASB(opt);
}


