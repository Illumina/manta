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

#include "options/CallOptionsDiploid.hh"


boost::program_options::options_description
getOptionsDescription(CallOptionsDiploid& opt)
{
    namespace po = boost::program_options;
    po::options_description desc("germline-variant-calling");
    desc.add_options()
    ("max-depth-factor", po::value(&opt.maxDepthFactor)->default_value(opt.maxDepthFactor),
     "Variants where the depth around the breakpoint is greater than this factor x the chromosomal mean will be filtered out")
    ("min-qual-score", po::value(&opt.minOutputAltScore)->default_value(opt.minOutputAltScore),
     "minimum QUAL score for variants included in the germline output vcf")
    ;

    return desc;
}

