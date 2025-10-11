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
