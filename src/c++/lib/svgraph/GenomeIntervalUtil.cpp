// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Manta
// Copyright (c) 2013-2014 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/sequencing/licenses/>
//

///
/// \author Chris Saunders
///

#include "svgraph/GenomeIntervalUtil.hh"



std::vector<unsigned>
intervalCompressor(
    std::vector<GenomeInterval>& intervals)
{
    std::vector<GenomeInterval> intervals2;

    const unsigned count(intervals.size());
    std::vector<bool> isTransfered(count,false);

    std::vector<unsigned> indexMap(count,0);

    for (unsigned headIndex(0); headIndex < count; ++headIndex)
    {
        if (isTransfered[headIndex]) continue;

        const unsigned headIndex2(intervals2.size());
        isTransfered[headIndex] = true;
        indexMap[headIndex] = headIndex2;
        intervals2.push_back(intervals[headIndex]);
        GenomeInterval& headInterval(intervals2.back());

        while (true)
        {
            bool isComplete(true);

            for (unsigned testIndex(headIndex+1); testIndex < count; ++testIndex)
            {
                if (isTransfered[testIndex]) continue;

                if (headInterval.isIntersect(intervals[testIndex]))
                {
                    isTransfered[testIndex] = true;
                    indexMap[testIndex] = headIndex2;
                    headInterval.range.merge_range(intervals[testIndex].range);
                    isComplete=false;
                    break;
                }
            }

            if (isComplete) break;
        }
    }

    intervals = intervals2;
    return indexMap;
}
