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

#include "alignment/AlignmentScores.hh"
#include "options/SmallAssemblerOptions.hh"


/// Options for the SV refiner step
///
/// Note that we have two categories of options for assembly an alignment,
/// one for small events, and one for large events
///
struct SVRefinerOptions
{
    /// match, mismatch, open score ratios taken from bwa defaults (but not extend!) :
    ///
    SVRefinerOptions() :
        smallSVAlignScores(2, -8, -12, 0, -1),
        spanningAlignScores(2, -8, -12, -1, -1),
        jumpScore(-25)
    {
        spanningAssembleOpt.minContigLength=75; ///< For breakend-spanning assemblies we require a larger contig than for small-variant assemblies
    }

    /// parameters for small SV assembly:
    AlignmentScores<int> smallSVAlignScores;
    SmallAssemblerOptions smallSVAssembleOpt;

    // parameters for large SV assembly:
    AlignmentScores<int> spanningAlignScores;
    int jumpScore;
    SmallAssemblerOptions spanningAssembleOpt;
};
