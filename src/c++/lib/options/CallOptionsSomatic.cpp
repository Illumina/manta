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

#include "options/CallOptionsSomatic.hh"


boost::program_options::options_description
getOptionsDescription(CallOptionsSomatic& opt)
{
    namespace po = boost::program_options;
    po::options_description desc("somatic-variant-calling");
    desc.add_options()
    ("somatic-max-depth-factor", po::value(&opt.maxDepthFactor)->default_value(opt.maxDepthFactor),
     "Variants where the normal-sample depth around the breakpoint is greater than this factor x the chromosomal mean will be filtered out")
    ("min-somatic-score", po::value(&opt.minOutputSomaticScore)->default_value(opt.minOutputSomaticScore),
     "minimum somatic quality score for variants to be included in the somatic output vcf")
    ("min-pass-somatic-score", po::value(&opt.minPassSomaticScore)->default_value(opt.minPassSomaticScore),
     "minimum somatic quality score below which variants are marked as filtered in the somatic output vcf")
    /*
        ("noise-sv-prior", po::value(&opt.noiseSVPrior)->default_value(opt.noiseSVPrior),
         "probability of a spurious SV observation shared in the tumor and normal samples")
    */
    ;

    return desc;
}
