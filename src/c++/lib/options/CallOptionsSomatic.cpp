//
// Manta - Structural Variant and Indel Caller
// Copyright (c) 2013-2025 Illumina, Inc.
//
// This program is licensed under the terms of the Polyform strict license
//
// ***As far as the law allows, the software comes as is, without
// any warranty or condition, and the licensor will not be liable
// to you for any damages arising out of these terms or the use
// or nature of the software, under any kind of legal claim.***
//
// You should have received a copy of the PolyForm Strict License 1.0.0
// along with this program.  If not, see <https://polyformproject.org/licenses/strict/1.0.0>.
//
//

#include "options/CallOptionsSomatic.hpp"

boost::program_options::options_description getOptionsDescription(CallOptionsSomatic& opt)
{
  namespace po = boost::program_options;
  po::options_description desc("somatic-variant-calling");
  // clang-format off
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
  // clang-format on

  return desc;
}
