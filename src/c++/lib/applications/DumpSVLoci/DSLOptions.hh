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

#pragma once

#include "manta/Program.hh"

#include <string>


struct DSLOptions
{
    DSLOptions() :
        isLocusIndex(false),
        locusIndex(0)
    {}

    bool isLocusIndex = false;
    unsigned locusIndex = 0;
    std::string graphFilename;
    std::string locusFilename;
    std::string region;
};


void
parseDSLOptions(const manta::Program& prog,
                int argc, char* argv[],
                DSLOptions& opt);
