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

#include "boost/program_options.hpp"


struct SomaticCallOptions
{

    SomaticCallOptions() :
        maxDepthFactor(3.0),
        minOutputSomaticScore(10)
    {}

    // breakpoints where the non-tumor depth is greater than the chromosome average x this factor are filtered out:
    float maxDepthFactor;

    unsigned minOutputSomaticScore;
};


boost::program_options::options_description
getOptionsDescription(SomaticCallOptions& opt);

