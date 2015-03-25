// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Manta
// Copyright (c) 2013-2015 Illumina, Inc.
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
    float germlineSVPrior = 1e-5;
    float somaticSVPrior = 1e-7;

    /// small/large values below reflect our expectation that there is more shared small event noise in small events
    ///
    /// expected shared tumor-normal sample noise rates for "small" SVs, ramp is from 3k->5k for small to large.
    float smallNoiseSVPrior = 1e-9;
    /// expected shared tumor-normal sample noise rates for "large" SVs
    float largeNoiseSVPrior = 1e-10;

    /// breakpoints where the non-tumor depth is greater than the chromosome average x this factor
    /// are filtered out:
    float maxDepthFactor = 3.0;
    std::string maxDepthFilterLabel = "MaxDepth";

    /// minimum somatic quality to print out a somatic variant
    unsigned minOutputSomaticScore = 10;
    /// minimum somatic quality which passes vcf filtration
    unsigned minPassSomaticScore = 30;
    std::string minSomaticScoreLabel = "MinSomaticScore";

    float maxMQ0Frac = 0.4;
    std::string maxMQ0FracLabel = "MaxMQ0Frac";
};


boost::program_options::options_description
getOptionsDescription(CallOptionsSomatic& opt);
