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

#pragma once

#include "manta/Program.hh"

#include <string>
#include <vector>



struct ESLOptions {

    ESLOptions() :
        minMapq(10),
        breakendEdgeTrimProb(0.2)
    {}

    /// report breakend regions with x prob regions removed from each edge
    unsigned minMapq;
    float breakendEdgeTrimProb;

    std::vector<std::string> alignmentFilename;
    std::string region;
    std::string statsFilename;
};


void
parseESLOptions(const manta::Program& prog,
                int argc, char* argv[],
                ESLOptions& opt);
