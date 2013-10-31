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
#include "blt_util/blt_types.hh"
#include "common/Exceptions.hh"

#include <cassert>
#include <cmath>
#include <iostream>


//#define DEBUG_SRA



std::ostream&
operator<<(std::ostream& os, const SRAlignmentInfo& info)
{
    os << "leftSize=" << info.leftSize << " rightSize=" << info.rightSize
       << " leftMismatches=" << info.leftMismatches << " rightMismatches=" << info.rightMismatches
       << " alignScore=" << info.alignScore
       << " isEvidence: " << info.isEvidence
       << " evidence: " << info.evidence
       << "\n";
    return os;
}



static
float
getLnLhood(
    const std::string& querySeq,
    const qscore_snp& qualConvert,
    const uint8_t* queryQual,
    const std::string::const_iterator& scanWindowBegin,
    const std::string::const_iterator& scanWindowEnd,
    const bool isBest,
    const float bestLnLhood)
{
    static const float ln_one_third(std::log(1./3.));

    const unsigned querySize(querySeq.size());

    float lnLhood(0);
    for (unsigned i(0); i<querySize; i++)
    {
        assert((scanWindowBegin+i) != scanWindowEnd);

        // put a lower-bound on quality values:
        const int baseQual(std::max(2,static_cast<int>(queryQual[i])));

        if ((querySeq[i] != *(scanWindowBegin+i)) || (querySeq[i] == 'N'))
        {
            if ((querySeq[i] == 'N') || (*(scanWindowBegin+i) == 'N'))
            {
                static const float lnRandomBase(std::log(0.25));
                lnLhood += lnRandomBase;
            }
            else
            {
                lnLhood += qualConvert.qphred_to_ln_error_prob(baseQual) + ln_one_third;
            }
        }
        else
        {
            lnLhood += qualConvert.qphred_to_ln_comp_error_prob(baseQual);
        }

        if (isBest && (lnLhood < bestLnLhood)) break;
    }

    return lnLhood;
}



static
void
calculateAlignScore(
    const std::string& querySeq,
    const std::string::const_iterator& scanWindowBegin,
    const std::string::const_iterator& scanWindowEnd,
    SRAlignmentInfo& alignment)
{
    const unsigned querySize = querySeq.size();
    alignment.leftMismatches=0;
    alignment.rightMismatches=0;

    for (unsigned i(0); i<querySize; i++)
    {
        assert((scanWindowBegin+i) != scanWindowEnd);

        if ((querySeq[i] != *(scanWindowBegin+i)) || (querySeq[i] == 'N'))
        {
            if (i<=alignment.leftSize)
            {
                alignment.leftMismatches += 1;
            }
            else
            {
                alignment.rightMismatches += 1;
            }
        }
    }

    alignment.alignScore = querySize - (alignment.leftMismatches+ alignment.rightMismatches);
}



static
void
setEvidence(
    SRAlignmentInfo& alignment)
{
    alignment.isEvidence = false;
    alignment.evidence = 0;

    const float size(static_cast<float>(alignment.leftSize+alignment.rightSize));

    //
    // filters for a read being counted as evidence
    //

    // adding new flank size threshold -- this might have to be changed based on sv size:
    static const unsigned minFlankSize(20);
    if (alignment.leftSize < minFlankSize) return;
    if (alignment.rightSize < minFlankSize) return;

    if ((alignment.leftMismatches/(float)alignment.leftSize) >= 0.25) return;
    if ((alignment.rightMismatches/(float)alignment.rightSize) >= 0.25) return;
    if ((alignment.alignScore/size) < 0.9) return;

    alignment.isEvidence = true;
    alignment.evidence = 2 * std::min(alignment.leftSize, alignment.rightSize) / (size);
}



void
splitReadAligner(
    const std::string& querySeq,
    const qscore_snp& qualConvert,
    const uint8_t* queryQual,
    const std::string& targetSeq,
    const unsigned targetBpBeginPos,
    SRAlignmentInfo& alignment)
{
    using namespace illumina::common;

    const unsigned querySize = querySeq.size();
    const unsigned targetSize = targetSeq.size();
    if (querySize >= targetSize)
    {
        std::ostringstream oss;
        oss << "ERROR: Unexpected split read alignment input."
            << " querySize: " << querySize << " targetSize: " << targetSize << '\n'
            << "\tquerySeq: " << querySeq << '\n'
            << "\ttargetSeq: " << targetSeq << '\n';
        BOOST_THROW_EXCEPTION(LogicException(oss.str()));

    }
    /// TODO: rerig things to use end_pos() as well

    // set the scanning start & end to make sure the candidate windows overlapping the breakpoint
    const unsigned scanStart(std::max(0, static_cast<pos_t>(targetBpBeginPos) - static_cast<pos_t>(querySize) + 2));
    const unsigned scanEnd(std::max(0, std::min((static_cast<pos_t>(targetBpBeginPos)-1), static_cast<pos_t>(targetSize - querySize)) ));

#ifdef DEBUG_SRA
    log_os << "query size = " << querySize << "target size = " << targetSize << "\n";
    log_os << "scan start = " << scanStart << " scan end = " << scanEnd << "\n";
#endif

    const std::string::const_iterator scanWindowBegin(targetSeq.begin());
    const std::string::const_iterator scanWindowEnd(targetSeq.end());

    // do one high-speed pass to find the optimal alignment (in terms of lhood), then compute all the goodies later:
    float bestLnLhood(0);
    unsigned bestPos(0);
    {
        bool isBest(false);
        for (unsigned i = scanStart; i<= scanEnd; i++)
        {
            const float lnLhood(getLnLhood(querySeq, qualConvert, queryQual,
                                           scanWindowBegin+i, scanWindowEnd, isBest, bestLnLhood));

            if ( (! isBest) || (lnLhood > bestLnLhood))
            {
                bestLnLhood = lnLhood;
                bestPos=i;
                isBest=true;
            }
        }
    }

    assert(bestPos <= (targetBpBeginPos+1));
    alignment.leftSize = static_cast<pos_t>(targetBpBeginPos+1) - bestPos;
    if (alignment.leftSize > querySize)
    {
        std::ostringstream oss;
        oss << "ERROR: Unexpected split read alignment outcome. "
            << " targetBeginPos: " << targetBpBeginPos << " bestPos: " << bestPos << " querySize: " << querySize << '\n'
            << "\tquerySeq: " << querySeq << '\n'
            << "\ttargetSeq: " << targetSeq << '\n';
        BOOST_THROW_EXCEPTION(LogicException(oss.str()));
    }
    alignment.rightSize= querySize - alignment.leftSize;
    alignment.alignLnLhood = bestLnLhood;
    alignment.alignPos = bestPos;

    calculateAlignScore(querySeq, scanWindowBegin+bestPos, scanWindowEnd, alignment);

    // filtering the alignment and set evidence
    setEvidence(alignment);

#ifdef DEBUG_SRA
    log_os << "final alignment\n" << alignment << "\n";
#endif
}
