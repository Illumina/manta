// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Manta
// Copyright (c) 2013-2014 Illumina, Inc.
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
    float indelPrior = 1e-5;

    /// breakpoints where the non-tumor depth is greater than the chromosome average x this factor are filtered out:
    float maxDepthFactor = 3.0;
    std::string maxDepthFilterLabel = "MaxDepth";

    /// minimum QUAL score to print out a diploid variant
    unsigned minOutputAltScore = 10;

    /// below this GQ value, the record is filtered in the diploid output VCF
    unsigned minPassGTScore = 20;
    std::string minGTFilterLabel = "MinGQ";

    // control filtration based on MQ0 fraction:
    float maxMQ0Frac = 0.4;
    std::string maxMQ0FracLabel = "MaxMQ0Frac";

    std::string rnaFilterLabel = "RNAFail";
};


boost::program_options::options_description
getOptionsDescription(CallOptionsDiploid& opt);
