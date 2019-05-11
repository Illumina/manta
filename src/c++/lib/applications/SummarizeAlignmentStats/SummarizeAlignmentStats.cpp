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

#include "SummarizeAlignmentStats.hpp"
#include "SASOptions.hpp"

#include "blt_util/log.hpp"
#include "common/OutStream.hpp"
#include "manta/ReadGroupStatsSet.hpp"

#include <iostream>

//#define OUTPUT_CDF

static void runSAS(const SASOptions& opt)
{
  static const float    quantLevel[] = {0.01f, 0.05f, 0.10f, 0.25f, 0.50f, 0.75f, 0.90f, 0.95f, 0.99f};
  static const unsigned quantLevelCount(sizeof(quantLevel) / sizeof(float));

  OutStream     outs(opt.outputFilename);
  std::ostream& report_os(outs.getStream());

  ReadGroupStatsSet rgss;
  rgss.load(opt.statsFilename.c_str());

  const unsigned groupCount(rgss.size());
  for (unsigned groupIndex(0); groupIndex < groupCount; ++groupIndex) {
    const ReadGroupStatsSet::KeyType& key(rgss.getKey(groupIndex));
#ifdef READ_GROUPS
    report_os << "bamFile:\t" << key.bamLabel << '\n';
    report_os << "readGroup:\t" << key.rgLabel << '\n';
#else
    report_os << "group:\t" << key.bamLabel << '\n';
#endif

    const ReadGroupStats& rgs(rgss.getStats(groupIndex));
    report_os << "fragment length observations:\t" << rgs.fragStats.totalObservations() << '\n';
    report_os << "fragment length quantiles:\n";
    for (unsigned quantLevelIndex(0); quantLevelIndex < quantLevelCount; ++quantLevelIndex) {
      report_os << quantLevel[quantLevelIndex] << '\t' << rgs.fragStats.quantile(quantLevel[quantLevelIndex])
                << '\n';
    }
    report_os << "Total sampled reads:\t" << rgs.readCounter.totalReadCount() << '\n';
    report_os << "Total sampled paired reads:\t" << rgs.readCounter.totalPairedReadCount() << '\n';
    report_os << "Total sampled unpaired reads:\t" << rgs.readCounter.totalUnpairedReadCount() << '\n';
    report_os << "Total sampled paired reads with low MAPQ:\t"
              << rgs.readCounter.totalPairedLowMapqReadCount() << '\n';
    report_os << "Total sampled high-confidence read pairs passing all filters:\t"
              << rgs.readCounter.totalHighConfidenceReadPairCount() << '\n';
    report_os << '\n';

#ifdef OUTPUT_CDF
    report_os << "fragment length cdf:\n";
    for (unsigned fragIndex(0); fragIndex <= 1000; ++fragIndex) {
      const unsigned fragLen(fragIndex * 10);
      report_os << fragLen << '\t' << rgs.fragStats.cdf(fragLen) << '\n';
    }
#endif
  }
}

void SummarizeAlignmentStats::runInternal(int argc, char* argv[]) const
{
  SASOptions opt;

  parseSASOptions(*this, argc, argv, opt);
  runSAS(opt);
}
