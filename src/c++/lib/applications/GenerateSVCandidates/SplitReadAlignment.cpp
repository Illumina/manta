// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Manta
// Copyright (c) 2013 Illumina, Inc.
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

#include "SplitReadAlignment.hh"
#include "blt_util/log.hh"
#include "blt_util/qscore.hh"

#include <cassert>
#include <cmath>


//#define DEBUG_SRA

unsigned
SplitReadAlignment::
calculateAlignScore(
    const std::string& querySeq,
    const uint8_t* queryQual,
    const std::string::const_iterator& scanWindowBegin,
    const std::string::const_iterator& scanWindowEnd)
{
    static const float ln_one_third(std::log(1./3.));

    const unsigned leftSize = _alignment.get_leftSize();
    const unsigned querySize = querySeq.size();
    unsigned leftMismatches(0);
    unsigned rightMismatches(0);

    float lnLhood(0);

    for (unsigned i(0); i<querySize; i++)
    {
        assert((scanWindowBegin+i) != scanWindowEnd);
        if (querySeq[i] != *(scanWindowBegin+i))
        {
            if (i<=leftSize)
            {
                leftMismatches += 1;
            }
            else
            {
                rightMismatches += 1;
            }
            lnLhood += qphred_to_ln_error_prob(static_cast<int>(queryQual[i])) + ln_one_third;
        }
        else
        {
            lnLhood += qphred_to_ln_comp_error_prob(static_cast<int>(queryQual[i]));
        }
    }

    const unsigned score = querySize - (leftMismatches+rightMismatches);
    _alignment.set_mismatches(leftMismatches, rightMismatches);
    _alignment.set_score(score);
    _alignment.set_lnLhood(lnLhood);

    return (score);
}

void
SplitReadAlignment::
align(const std::string& querySeq,
      const uint8_t* queryQual,
      const std::string& targetSeq,
      const unsigned bpOffset)
{
    SRAlignmentInfo bestAlignInfo;

    const unsigned querySize = querySeq.size();
    const unsigned targetSize = targetSeq.size();
    assert(querySize < targetSize);

    // set the scanning start & end to make sure the candidate windows overlapping the breakpoint
    const unsigned scanStart = ((bpOffset+2) <= querySize)? 0 : (bpOffset - querySize + 2);
    const unsigned scanEnd = (bpOffset == 0)? 0 : std::min((bpOffset-1), (targetSize - querySize));

#ifdef DEBUG_SRA
    log_os << "query size = " << querySize << "target size = " << targetSize << "\n";
    log_os << "scan start = " << scanStart << " scan end = " << scanEnd << "\n";
#endif

    const std::string::const_iterator scanWindowEnd(targetSeq.end());

    for (unsigned i = scanStart; i<= scanEnd; i++)
    {
        const unsigned leftSize = bpOffset - i + 1;
        const unsigned rightSize= querySize - leftSize;
        _alignment.set_sizes(leftSize, rightSize);

        const std::string::const_iterator scanWindowStart = targetSeq.begin()+i;
        const unsigned score = calculateAlignScore(querySeq, queryQual, scanWindowStart, scanWindowEnd);

        if (score > bestAlignInfo.get_alignScore())
        {
            bestAlignInfo.update_alignInfo(_alignment);
        }
    }

    // update the alignment with the best one
    _alignment.update_alignInfo(bestAlignInfo);
    // filtering the alignment and set evidence
    set_evidence();

#ifdef DEBUG_SRA
    log_os << "final alignment\n" << _alignment << "\n";
#endif
}


void
SplitReadAlignment::
set_evidence()
{
    _hasEvidence = false;
    const unsigned leftSize = _alignment.get_leftSize();
    const unsigned rightSize = _alignment.get_rightSize();
    const unsigned bestScore = leftSize + rightSize;

    // filters for a read being counted as evidence
    if (((_alignment.get_leftMismatches()/(float)leftSize) < 0.25) &&
        ((_alignment.get_rightMismatches()/(float)rightSize) < 0.25) &&
        ((_alignment.get_alignScore()/(float)bestScore) >= 0.9))
    {
        _hasEvidence = true;
        _evidence = 2 * std::min(leftSize, rightSize) / (float)(leftSize + rightSize);
    }
}

std::ostream&
operator<<(std::ostream& os, const SRAlignmentInfo& info)
{
    os << "leftSize=" << info.get_leftSize() << " rightSize=" << info.get_rightSize()
       << " leftMismatches=" << info.get_leftMismatches() << " rightMismatches=" << info.get_rightMismatches()
       << " alignScore=" << info.get_alignScore() << "\n";
    return os;
}

std::ostream&
operator<<(std::ostream& os, const SplitReadAlignment& srAlign)
{
    os << "has_evidence=" << srAlign.has_evidence() << " evidence=" << srAlign.get_evidence() << "\n";
    os << "alignment:\n" << srAlign.get_alignment();
    return os;
}
