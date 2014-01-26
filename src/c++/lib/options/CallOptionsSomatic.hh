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


struct CallOptionsSomatic
{
    CallOptionsSomatic() :
        germlineSVPrior(1e-5),
        somaticSVPrior(1e-7),
        maxDepthFactor(3.0),
        maxDepthFilterLabel("MaxDepth"),
        minOutputSomaticScore(10),
        minPassSomaticScore(30),
        minSomaticScoreLabel("MinSomaticScore"),
        maxMQ0Frac(0.4),
        maxMQ0FracLabel("MaxMQ0Frac")

    {}

    float germlineSVPrior;
    float somaticSVPrior;

    // breakpoints where the non-tumor depth is greater than the chromosome average x this factor are filtered out:
    float maxDepthFactor;
    std::string maxDepthFilterLabel;

    unsigned minOutputSomaticScore; ///< minimum somatic quality to print out a somatic variant
    unsigned minPassSomaticScore; ///< minimum somatic quality which passes vcf filtration
    std::string minSomaticScoreLabel;

    float maxMQ0Frac;
    std::string maxMQ0FracLabel;
};


boost::program_options::options_description
getOptionsDescription(CallOptionsSomatic& opt);
