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

#include "GetChromDepth.hpp"
#include "ChromDepthOptions.hpp"
#include "ReadChromDepthUtil.hpp"

#include "blt_util/log.hpp"
#include "common/OutStream.hpp"

#include <cstdlib>

#include <iomanip>
#include <iostream>

static void getChromDepth(const ChromDepthOptions& opt)
{
  // check that we have write permission on the output file early:
  {
    OutStream outs(opt.outputFilename);
  }

  std::vector<double> chromDepth;
  for (const std::string& chromName : opt.chromNames) {
    chromDepth.push_back(
        readChromDepthFromAlignment(opt.referenceFilename, opt.alignmentFilename, chromName));
  }

  OutStream     outs(opt.outputFilename);
  std::ostream& os(outs.getStream());

  const unsigned chromCount(opt.chromNames.size());
  for (unsigned chromIndex(0); chromIndex < chromCount; ++chromIndex) {
    os << opt.chromNames[chromIndex] << "\t" << std::fixed << std::setprecision(2) << chromDepth[chromIndex]
       << "\n";
  }
}

void GetChromDepth::runInternal(int argc, char* argv[]) const
{
  ChromDepthOptions opt;

  parseChromDepthOptions(*this, argc, argv, opt);
  getChromDepth(opt);
}
