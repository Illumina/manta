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
#include "options/ReadScannerOptions.hh"

#include <string>
#include <vector>



struct GSCOptions
{
    GSCOptions() :
        binCount(1),
        binIndex(0)
    {}

    ReadScannerOptions scanOpt;

    std::vector<std::string> alignmentFilename;
    std::vector<bool> isAlignmentTumor;
    std::string graphFilename;
    std::string referenceFilename;
    std::string statsFilename;

    std::string candidateOutputFilename;
    //std::string germlineOutputFilename;
    std::string somaticOutputFilename;

    unsigned binCount;
    unsigned binIndex;
};


void
parseGSCOptions(const manta::Program& prog,
                int argc, char* argv[],
                GSCOptions& opt);
