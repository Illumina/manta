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

#include "EstimateSVLoci.hpp"
#include "EstimateSVLociRunner.hpp"

#include "common/OutStream.hpp"

static void runEstimateSVLoci(const ESLOptions& opt)
{
  {
    // early test that we have permission to write to output file
    OutStream outs(opt.outputFilename);
  }

  EstimateSVLociRunner eslRunner(opt);
  for (const auto& region : opt.regions) {
    eslRunner.estimateSVLociForSingleRegion(region);
  }

  eslRunner.getLocusSet().save(opt.outputFilename.c_str());
}

void EstimateSVLoci::runInternal(int argc, char* argv[]) const
{
  ESLOptions opt;

  parseESLOptions(*this, argc, argv, opt);
  runEstimateSVLoci(opt);
}
