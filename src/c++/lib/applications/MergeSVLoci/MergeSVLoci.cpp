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

#include "MergeSVLoci.hpp"
#include "MSLOptions.hpp"

#include "blt_util/log.hpp"
#include "common/OutStream.hpp"
#include "svgraph/SVLocusSet.hpp"

static void runMSL(const MSLOptions& opt)
{
  TimeTracker timer;
  timer.resume();
  {
    // early test that we have permission to write to output file
    OutStream outs(opt.outputFilename);
  }

  const unsigned graphFileCount(opt.graphFilename.size());
  // This should already be enforced by the arg parsing interface:
  assert(graphFileCount > 0);

  if (opt.isVerbose) {
    log_os << "INFO: Initializing from file: '" << opt.graphFilename[0] << "'\n";
  }

  SVLocusSet mergedSet(opt.graphFilename[0].c_str());

  if (opt.isVerbose) {
    log_os << "INFO: Finished initializing from file: '" << opt.graphFilename[0] << "'\n";
  }

  for (unsigned graphFileIndex(1); graphFileIndex < graphFileCount; ++graphFileIndex) {
    const std::string& graphFile(opt.graphFilename[graphFileIndex]);

    if (opt.isVerbose) {
      log_os << "INFO: Merging file: '" << graphFile << "'\n";
    }

    const SVLocusSet inputSet(graphFile.c_str());
    mergedSet.merge(inputSet);

    if (opt.isVerbose) {
      log_os << "INFO: Finished merging file: '" << graphFile << "'\n";
    }
  }

  mergedSet.finalize();
  if (opt.isVerbose) {
    log_os << "INFO: Finished cleaning merged graph.\n";
  }
  timer.stop();
  mergedSet.addMergeTime(timer.getTimes());
  mergedSet.save(opt.outputFilename.c_str());
}

void MergeSVLoci::runInternal(int argc, char* argv[]) const
{
  MSLOptions opt;

  parseMSLOptions(*this, argc, argv, opt);
  runMSL(opt);
}
