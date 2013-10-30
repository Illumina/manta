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

#include "boost/program_options.hpp"


struct CallOptionsDiploid
{

    CallOptionsDiploid() :
        indelPrior(1e-5),
        maxDepthFactor(3.0),
        maxDepthFilterLabel("MaxDepth"),
        minOutputAltScore(10),
        minGTScoreFilter(20),
        minGTFilterLabel("MinGQ")
    {}

    float indelPrior;

    // breakpoints where the non-tumor depth is greater than the chromosome average x this factor are filtered out:
    float maxDepthFactor;
    std::string maxDepthFilterLabel;

    unsigned minOutputAltScore; ///< minimum QUAL score to print out a diploid variant

    unsigned minGTScoreFilter; ///< below this GQ value, the record is filtered in the diploid output VCF
    std::string minGTFilterLabel;
};


boost::program_options::options_description
getOptionsDescription(CallOptionsDiploid& opt);
