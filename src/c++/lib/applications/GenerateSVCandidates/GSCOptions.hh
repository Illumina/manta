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

#include "EdgeOptions.hh"
#include "manta/Program.hh"
#include "options/AlignmentFileOptions.hh"
#include "options/ReadScannerOptions.hh"
#include "options/SomaticCallOptions.hh"

#include <string>
#include <vector>


struct GSCOptions
{
    GSCOptions() :
        isSkipAssembly(false)
    {}

    AlignmentFileOptions alignFileOpt;
    EdgeOptions edgeOpt;
    ReadScannerOptions scanOpt;
    SomaticCallOptions somaticOpt;

    std::string graphFilename;
    std::string referenceFilename;
    std::string statsFilename;
    std::string chromDepthFilename;

    std::string candidateOutputFilename;
    std::string somaticOutputFilename;

    bool isSkipAssembly; ///< if true, skip assembly and run a low-resolution, breakdancer-like subset of the workflow
};


void
parseGSCOptions(const manta::Program& prog,
                int argc, char* argv[],
                GSCOptions& opt);
