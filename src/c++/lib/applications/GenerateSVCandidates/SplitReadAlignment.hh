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
/// \author Xiaoyu Chen
///

#pragma once

#include "blt_util/known_pos_range2.hh"
#include "blt_util/qscore_snp.hh"

#include <stdint.h>

#include <string>
#include <iosfwd>


struct SRAlignmentInfo
{
    SRAlignmentInfo():
        alignPos(0),
        leftSize(0),
        homSize(0),
        rightSize(0),
        leftMismatches(0),
        homMismatches(0),
        rightMismatches(0),
        alignScore(0),
        alignLnLhood(0),
        isEvidence(false),
        isTier2Evidence(false),
        evidence(0)
    {}

    unsigned
    leftSplitSize(
        const bool isSplitLeftOfHomologyRange) const
    {
        return (leftSize + (isSplitLeftOfHomologyRange ? 0 : homSize));
    }

    unsigned
    rightSplitSize(
        const bool isSplitLeftOfHomologyRange) const
    {
        return (rightSize + (isSplitLeftOfHomologyRange ? homSize : 0));
    }

    unsigned
    leftSplitMismatches(
        const bool isSplitLeftOfHomologyRange) const
    {
        return (leftMismatches + (isSplitLeftOfHomologyRange ? 0 : homMismatches));
    }

    unsigned
    rightSplitMismatches(
        const bool isSplitLeftOfHomologyRange) const
    {
        return (rightMismatches + (isSplitLeftOfHomologyRange ? homMismatches : 0));
    }

    unsigned alignPos;
    unsigned leftSize;
    unsigned homSize;
    unsigned rightSize;
    unsigned leftMismatches;
    unsigned homMismatches;
    unsigned rightMismatches;
    unsigned alignScore;
    float alignLnLhood;

    bool isEvidence;
    bool isTier2Evidence;
    float evidence;
};

std::ostream&
operator<<(std::ostream& os, const SRAlignmentInfo& info);


///
/// \param[in] flankScoreSize the number of bases to score past the end of homology range
///
/// \param[in] targetBpOffsetRange this is the range of the breakend (accounting for homology) in targetSeq coordinates
///
/// \param[in] isSplitLeftOfHomologyRange for each split problem we need to split to either the left or the right end of the
///            breakend homology range
///
/// TODO: need to add a query subset/length limit, so that as the query size goes up (ie. 2 x 400) we still consistently
///       detect split read support without having to add more and more reference to the targetSeq
///
void
splitReadAligner(
    const unsigned flankScoreSize,
    const std::string& querySeq,
    const qscore_snp& qualConvert,
    const uint8_t* queryQual,
    const std::string& targetSeq,
    const known_pos_range2& targetBpOffsetRange,
    const bool isSplitLeftOfHomologyRange,
    SRAlignmentInfo& alignment);
