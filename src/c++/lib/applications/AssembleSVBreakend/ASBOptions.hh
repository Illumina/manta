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
// <https://github.com/sequencing/licenses/>
//

///
/// \author Ole Schulz-Trieglaff
///

#pragma once

#include "manta/Program.hh"
#include "manta/SVCandidate.hh"

#include "applications/GenerateSVCandidates/GSCOptions.hh"

#include <string>



struct ASBOptions : public GSCOptions
{

    ASBOptions() :
        GSCOptions(),
        breakend1(""),
        breakend2(""),
        contigOutfile("")
    {}

    std::string breakend1;
    std::string breakend2;
    std::string contigOutfile;
    std::string referenceFilename;

};


void
parseASBOptions(const manta::Program& prog,
                int argc, char* argv[],
                ASBOptions& opt);
