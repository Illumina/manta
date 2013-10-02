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

#include "splitReadAlignment.hh"

unsigned
splitReadAlignment::
calculateAlignScore(
    const std::string& querySeq,
    const std::string& scanWindow)
{
    const unsigned leftSize = _alignment.get_leftSize();
    const unsigned querySize = querySeq.size();
    unsigned leftMismatches(0);
    unsigned rightMismatches(0);

    for (unsigned i = 0; i<= leftSize; i++)
        if (querySeq[i] != scanWindow[i]) leftMismatches += 1;

    for (unsigned i = leftSize+1; i<querySize; i++)
        if (querySeq[i] != scanWindow[i]) rightMismatches += 1;

    const unsigned score = querySize - (leftMismatches+rightMismatches);
    _alignment.set_mismatches(leftMismatches, rightMismatches);
    _alignment.set_score(score);

    return (score);
}

void
splitReadAlignment::
align(const std::string& querySeq,
      const std::string& targetSeq,
      const unsigned bpOffset)
{
    SRAlignmentInfo bestAlignInfo;

    const unsigned querySize = querySeq.size();
    const unsigned targetSize = targetSeq.size();
    // set the scanning start & end to make sure the candidate windows overlapping the breakpoint
    const unsigned scanStart = std::max((unsigned)0, (bpOffset - querySize + 2));
    const unsigned scanEnd = std::min((bpOffset-1), (targetSize - querySize));

    for (unsigned i = scanStart; i<= scanEnd; i++)
    {
        const unsigned leftSize = bpOffset - i + 1;
        const unsigned rightSize= querySize - leftSize;
        _alignment.set_sizes(leftSize, rightSize);

        const std::string scanWindow = targetSeq.substr(i, querySize);
        const unsigned score = calculateAlignScore(querySeq, scanWindow);

        if (score > bestAlignInfo.get_alignScore())
        {
            bestAlignInfo.update_alignInfo(_alignment);
        }
    }

    // update the alignment with the best one
    _alignment.update_alignInfo(bestAlignInfo);
    // filtering the alignment and set evidence
    set_evidence();
}


void
splitReadAlignment::
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
operator<<(std::ostream& os, const splitReadAlignment& srAlign)
{
    os << "has_evidence=" << srAlign.has_evidence() << " evidence=" << srAlign.get_evidence() << "\n";
    os << "alignment:\n" << srAlign.get_alignment();
    return os;
}
