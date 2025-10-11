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

#include "options/ReadScannerOptionsParser.hpp"

boost::program_options::options_description getOptionsDescription(ReadScannerOptions& opt)
{
  namespace po = boost::program_options;
  po::options_description desc("read-scanner");
  // clang-format off
  desc.add_options()
  ("min-candidate-sv-size", po::value(&opt.minCandidateVariantSize)->default_value(opt.minCandidateVariantSize),
   "Indels below this size will not be discovered or reported as candidates")
  ("min-mapq", po::value(&opt.minMapq)->default_value(opt.minMapq),
   "Reads with MAPQ less than this value will be ignored")
  ("edge-prob", po::value(&opt.breakendEdgeQuantileProb)->default_value(opt.breakendEdgeQuantileProb),
   "Breakend range associated with each read will trimmed to expected fragment quantile range [p,(1-p)], p: edge-prob")
  ("use-overlapping-pair", po::value(&opt.useOverlapPairEvidence)->zero_tokens(),
   "Consider an overlapping read pair as evidence")
  ("ignore-anom-proper-pair", po::value(&opt.isIgnoreAnomProperPair)->zero_tokens(),
   "Disregard anomalous fragment sizes if the BAM record has the proper pair bit set. "
   "This flag is typically set for RNA-SEQ analysis where the proper-pair bit is used to indicate an intron-spanning read pair.")
  ;
  // clang-format on

  return desc;
}

bool parseOptions(
    const boost::program_options::variables_map& /*vm*/, ReadScannerOptions& opt, std::string& errorMsg)
{
  errorMsg.clear();
  if ((opt.breakendEdgeQuantileProb <= 0) || (opt.breakendEdgeQuantileProb >= 1.0)) {
    errorMsg = "edge-prob argument is restricted to (0,1)";
  }

  return (!errorMsg.empty());
}
