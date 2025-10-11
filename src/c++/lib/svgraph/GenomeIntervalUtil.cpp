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

#include "GenomeIntervalUtil.hpp"
#include "htsapi/bam_header_info.hpp"
#include "htsapi/bam_header_util.hpp"

std::vector<unsigned> intervalCompressor(std::vector<GenomeInterval>& intervals)
{
  std::vector<GenomeInterval> intervals2;

  const unsigned    count(intervals.size());
  std::vector<bool> isTransfered(count, false);

  std::vector<unsigned> indexMap(count, 0);

  for (unsigned headIndex(0); headIndex < count; ++headIndex) {
    if (isTransfered[headIndex]) continue;

    const unsigned headIndex2(intervals2.size());
    isTransfered[headIndex] = true;
    indexMap[headIndex]     = headIndex2;
    intervals2.push_back(intervals[headIndex]);
    GenomeInterval& headInterval(intervals2.back());

    while (true) {
      bool isComplete(true);

      for (unsigned testIndex(headIndex + 1); testIndex < count; ++testIndex) {
        if (isTransfered[testIndex]) continue;

        if (headInterval.isIntersect(intervals[testIndex])) {
          isTransfered[testIndex] = true;
          indexMap[testIndex]     = headIndex2;
          headInterval.range.merge_range(intervals[testIndex].range);
          isComplete = false;
          break;
        }
      }

      if (isComplete) break;
    }
  }

  intervals = intervals2;
  return indexMap;
}

GenomeInterval convertSamtoolsRegionToGenomeInterval(
    const bam_header_info& bamHeader, const std::string& region)
{
  int32_t tid(0), beginPos(0), endPos(0);
  parse_bam_region(bamHeader, region.c_str(), tid, beginPos, endPos);
  return GenomeInterval(tid, beginPos, endPos);
}
