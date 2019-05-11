//
// Manta - Structural Variant and Indel Caller
// Copyright (c) 2013-2019 Illumina, Inc.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
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
