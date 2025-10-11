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

#include "MergeAlignmentStatsOptions.hpp"

#include "blt_util/log.hpp"
#include "common/OutStream.hpp"
#include "manta/ReadGroupStatsUtil.hpp"

#include <cstdlib>

#include <iostream>
#include "MergeAlignmentStats.hpp"

static void mergeAlignmentStats(const MergeAlignmentStatsOptions& opt)
{
  if (opt.statsFiles.empty()) {
    log_os << "ERROR: No input files specified.\n";
    exit(EXIT_FAILURE);
  }

  ReadGroupStatsSet all_rstats;
  for (const std::string& file : opt.statsFiles) {
    ReadGroupStatsSet rstats;
    rstats.load(file.c_str());
    all_rstats.merge(rstats);
  }

  all_rstats.save(opt.outputFilename.c_str());
}

void MergeAlignmentStats::runInternal(int argc, char* argv[]) const
{
  MergeAlignmentStatsOptions opt;

  parseMergeAlignmentStatsOptions(*this, argc, argv, opt);
  mergeAlignmentStats(opt);
}
