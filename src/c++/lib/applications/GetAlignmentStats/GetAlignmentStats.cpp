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

#include "GetAlignmentStats.hpp"

#include "AlignmentStatsOptions.hpp"

#include "blt_util/log.hpp"
#include "common/OutStream.hpp"
#include "manta/ReadGroupStatsUtil.hpp"

#include <cstdlib>

#include <iostream>

static void runAlignmentStats(const AlignmentStatsOptions& opt)
{
  // calculate fragment size statistics for all read groups in all bams

  // instantiate early to test for filename/permissions problems
  if (opt.alignFileOpt.alignmentFilenames.empty()) {
    log_os << "ERROR: No input files specified.\n";
    exit(EXIT_FAILURE);
  }

  ReadGroupStatsSet rstats;

  for (const std::string& alignmentFilename : opt.alignFileOpt.alignmentFilenames) {
    extractReadGroupStatsFromAlignmentFile(
        opt.referenceFilename, alignmentFilename, opt.defaultStatsFilename, rstats);
  }

  rstats.save(opt.outputFilename.c_str());
}

void GetAlignmentStats::runInternal(int argc, char* argv[]) const
{
  AlignmentStatsOptions opt;

  parseAlignmentStatsOptions(*this, argc, argv, opt);
  runAlignmentStats(opt);
}
