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

/// \file
/// \author Chris Saunders
///

#include "DumpSVLoci.hpp"
#include "DSLOptions.hpp"

#include "svgraph/GenomeIntervalUtil.hpp"
#include "svgraph/SVLocusSet.hpp"

#include "blt_util/thirdparty_push.h"

#include "boost/archive/binary_oarchive.hpp"

#include "blt_util/thirdparty_pop.h"

#include <fstream>
#include <iostream>

static void runDSL(const DSLOptions& opt)
{
  SVLocusSet set(opt.graphFilename.c_str());

  const SVLocusSet& cset(set);

  std::ostream& os(std::cout);

  // add this handy map of chromosome id to chromosome label at the start of all output types:
  os << cset.getBamHeader() << "\n";

  if (!opt.region.empty()) {
    set.dumpRegion(os, convertSamtoolsRegionToGenomeInterval(cset.getBamHeader(), opt.region));
  } else if (opt.isLocusIndex) {
    const SVLocus& locus(cset.getLocus(opt.locusIndex));
    if (opt.locusFilename.empty()) {
      os << locus;
    } else {
      std::ofstream                   ofs(opt.locusFilename.c_str(), std::ios::binary);
      boost::archive::binary_oarchive oa(ofs);
      oa << locus;
    }
  } else {
    cset.dump(os);
  }
}

void DumpSVLoci::runInternal(int argc, char* argv[]) const
{
  DSLOptions opt;

  parseDSLOptions(*this, argc, argv, opt);
  runDSL(opt);
}
