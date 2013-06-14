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


#pragma once

#include "manta/Program.hh"

#include <string>
#include <vector>



struct AlignmentStatsOptions {

    std::vector<std::string> alignmentFilename;
};


void
parseAlignmentStatsOptions(const manta::Program& prog,
                           int argc, char* argv[],
                           AlignmentStatsOptions& opt);
