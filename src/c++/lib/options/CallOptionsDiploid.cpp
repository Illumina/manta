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

#include "options/CallOptionsDiploid.hpp"

boost::program_options::options_description getOptionsDescription(CallOptionsDiploid& opt)
{
  namespace po = boost::program_options;
  po::options_description desc("germline-variant-calling");
  // clang-format off
  desc.add_options()
  ("diploid-max-depth-factor", po::value(&opt.maxDepthFactor)->default_value(opt.maxDepthFactor),
   "Variants where the depth around the breakpoint is greater than this factor x the chromosomal mean will be filtered out")
  ("min-qual-score", po::value(&opt.minOutputAltScore)->default_value(opt.minOutputAltScore),
   "minimum QUAL score for variants included in the germline output vcf")
  ("min-pass-qual-score", po::value(&opt.minPassAltScore)->default_value(opt.minPassAltScore),
   "minimum QUAL score for variants to PASS in germline output vcf")
  ("min-pass-gt-score", po::value(&opt.minPassGTScore)->default_value(opt.minPassGTScore),
   "minimum genotype quality score below which samples are filtered for a variant in the germline output vcf")
  ;
  // clang-format on

  return desc;
}
