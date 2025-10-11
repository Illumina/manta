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

#include "options/CallOptionsTumor.hpp"

boost::program_options::options_description getOptionsDescription(CallOptionsTumor& opt)
{
  namespace po = boost::program_options;
  po::options_description desc("tumor-only-variant-calling");
  // clang-format off
  desc.add_options()
  ("tumor-max-depth-factor", po::value(&opt.maxDepthFactor)->default_value(opt.maxDepthFactor),
   "Variants where the tumor-sample depth around the breakpoint is greater than this factor x the chromosomal mean will be filtered out")
  ;
  // clang-format on

  return desc;
}
