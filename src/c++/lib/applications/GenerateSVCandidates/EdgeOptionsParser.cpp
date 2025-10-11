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

/// \file
/// \author Chris Saunders
///

#include "EdgeOptionsParser.hpp"

#include "blt_util/parse_util.hpp"
#include "blt_util/string_util.hpp"

namespace {
const char locusIndexKey[] = "locus-index";
}

boost::program_options::options_description getOptionsDescription(EdgeOptions& opt)
{
  namespace po = boost::program_options;

  po::options_description optdesc("edge-selection");
  // clang-format off
  optdesc.add_options()
  ("bin-count", po::value(&opt.binCount)->default_value(opt.binCount),
   "Specify how many bins the SV candidate problem should be divided into, where bin-index can be used to specify which bin to solve")
  ("bin-index", po::value(&opt.binIndex)->default_value(opt.binIndex),
   "Specify which bin to solve when the SV candidate problem is subdivided into bins. Value must bin in [0,bin-count)")
  ("max-edge-count", po::value(&opt.graphNodeMaxEdgeCount)->default_value(opt.graphNodeMaxEdgeCount),
   "Specify the maximum number of edge count. If both nodes of an edge have an edge count higher than this, the edge is skipped for evaluation.")
  (locusIndexKey, po::value<std::string>(),
   "Instead of solving for all SV candidates in a bin, solve for candidates of a particular locus or edge."
   " If this argument is specified then bin-index is ignored."
   " Argument can be one of { locusIndex , locusIndex:nodeIndex , locusIndex:nodeIndex:nodeIndex },"
   " which will run an entire locus, all edges connected to one node in a locus or a single edge, respectively.")
  ;
  // clang-format on

  return optdesc;
}

bool parseOptions(const boost::program_options::variables_map& vm, EdgeOptions& opt, std::string& errorMsg)
{
  errorMsg.clear();

  if (vm.count(locusIndexKey)) {
    using namespace illumina::blt_util;

    const std::string& locusString(boost::any_cast<std::string>(vm[locusIndexKey].value()));

    std::vector<std::string> indices;
    split_string(locusString, ':', indices);
    if (indices.size() > 3) {
      errorMsg = "locus-index argument can have no more than 3 colon separated segments";
    }

    assert(!indices.empty());

    opt.isLocusIndex = true;
    {
      LocusEdgeOptions& lopt(opt.locusOpt);
      lopt.locusIndex = parse_unsigned_str(indices[0]);
      if (indices.size() > 1) {
        lopt.isNodeIndex1 = true;
        lopt.nodeIndex1   = parse_unsigned_str(indices[1]);
        if (indices.size() > 2) {
          lopt.isNodeIndex2 = true;
          lopt.nodeIndex2   = parse_unsigned_str(indices[2]);
        }
      }
    }
  }

  if (errorMsg.empty()) {
    if (opt.binCount < 1) {
      errorMsg = "bin-count must be 1 or greater";
    } else if (opt.binIndex >= opt.binCount) {
      errorMsg = "bin-index must be in range [0,bin-count)";
    }
  }

  return (!errorMsg.empty());
}
