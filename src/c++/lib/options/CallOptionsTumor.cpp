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
/// \author Xiaoyu Chen
///

#include "options/CallOptionsTumor.hh"


boost::program_options::options_description
getOptionsDescription(CallOptionsTumor& opt)
{
    namespace po = boost::program_options;
    po::options_description desc("tumor-only-variant-calling");
    desc.add_options()
    ("tumor-max-depth-factor", po::value(&opt.maxDepthFactor)->default_value(opt.maxDepthFactor),
      "Variants where the tumor-sample depth around the breakpoint is greater than this factor x the chromosomal mean will be filtered out")
    ;

    return desc;
}





