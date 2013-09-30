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
/// \author Chris Saunders
///

#pragma once

#include "EdgeOptions.hh"
#include "manta/Program.hh"
#include "options/AlignmentFileOptions.hh"
#include "options/ReadScannerOptions.hh"
#include "options/SomaticCallOptions.hh"
#include "options/SVRefinerOptions.hh"

#include <string>
#include <vector>


struct GSCOptions
{
    GSCOptions() :
        isSkipAssembly(false),
        minScoredVariantSize(51)
    {}

    AlignmentFileOptions alignFileOpt;
    EdgeOptions edgeOpt;
    ReadScannerOptions scanOpt;
    SVRefinerOptions refineOpt;
    SomaticCallOptions somaticOpt;

    std::string graphFilename;
    std::string referenceFilename;
    std::string statsFilename;
    std::string chromDepthFilename;

    std::string candidateOutputFilename;
    std::string somaticOutputFilename;

    bool isSkipAssembly; ///< if true, skip assembly and run a low-resolution, breakdancer-like subset of the workflow

    unsigned minScoredVariantSize; ///< min size for scoring and scored output following candidate generation
};


void
parseGSCOptions(const manta::Program& prog,
                int argc, char* argv[],
                GSCOptions& opt);
