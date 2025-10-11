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

#include "options/SVLocusSetOptionsParser.hpp"

boost::program_options::options_description getOptionsDescription(SVLocusSetOptions& opt)
{
  namespace po = boost::program_options;
  po::options_description desc("sv-locus-graph");
  // clang-format off
  desc.add_options()
  ("min-edge-observations", po::value(&opt.minMergeEdgeObservations)->default_value(opt.minMergeEdgeObservations),
   "Minimum number of supporting observations required to retain a graph edge")
  ;
  // clang-format on

  return desc;
}

bool parseOptions(
    const boost::program_options::variables_map& /*vm*/, SVLocusSetOptions& /*opt*/, std::string& errorMsg)
{
  errorMsg.clear();
#if 0
    if ((opt.breakendEdgeQuantileProb <= 0) || (opt.breakendEdgeQuantileProb >= 1.0))
    {
        errorMsg="edge-prob argument is restricted to (0,1)";
    }
#endif
  return (!errorMsg.empty());
}
