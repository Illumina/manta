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

#include "manta/Program.hh"
#include "manta/SVLocusScanner.hh"
#include "options/AlignmentFileOptions.hh"
#include "options/ReadScannerOptions.hh"
#include "options/SVLocusSetOptions.hh"


struct ESLOptions
{
    ESLOptions(const unsigned minMergeEdgeObservations = 3) :
        graphOpt(minMergeEdgeObservations * SVObservationWeights::observation),   // initialize noise edge filtration parameters
        maxDepthFactor(12.0)
    {}

    AlignmentFileOptions alignFileOpt;
    ReadScannerOptions scanOpt;
    SVLocusSetOptions graphOpt;

    std::string referenceFilename;
    std::string outputFilename;
    std::string region;
    std::string statsFilename;
    std::string chromDepthFilename;
    std::string truthVcfFilename;

    float maxDepthFactor; ///< the maximum depth at which information is entered into the graph when avg chrom depth is given
};


void
parseESLOptions(const manta::Program& prog,
                int argc, char* argv[],
                ESLOptions& opt);
