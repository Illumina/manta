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
#include "blt_util/blt_types.hh"
#include "blt_util/log.hh"
#include "blt_util/seq_printer.hh"
#include "common/Exceptions.hh"

#include <cassert>
#include <cmath>
#include <iostream>


//#define DEBUG_SRA



std::ostream&
operator<<(std::ostream& os, const SRAlignmentInfo& info)
{
    os << "leftSize=" << info.leftSize << " homSize=" << info.homSize << " rightSize=" << info.rightSize
       << " leftMismatches=" << info.leftMismatches << " homMismatches=" << info.homMismatches << " rightMismatches=" << info.rightMismatches
       << " alignScore=" << info.alignScore
       << " isEvidence: " << info.isEvidence
       << " evidence: " << info.evidence
       << " alignLnLhood: "  << info.alignLnLhood
       << '\n';
    return os;
}



/// \params[out] return the LnLhood expected from a perfect match to the reference
static
float
getLnLhood(
    const std::string& querySeq,
    const qscore_snp& qualConvert,
    const uint8_t* queryQual,
    const std::string& targetSeq,
    const pos_t targetStartOffset,
    const known_pos_range2& scoreRange,
    const bool isBest,
    const float bestLnLhood)
{
    static const float ln_one_third(std::log(1./3.));

    const unsigned querySize(querySeq.size());

    assert((targetStartOffset+querySize) <= targetSeq.size());

    float lnLhood(0);
    for (unsigned i(0); i<querySize; i++)
    {
        // put a lower-bound on quality values:
        const int baseQual(std::max(2,static_cast<int>(queryQual[i])));

        if ((targetStartOffset+static_cast<pos_t>(i)) > scoreRange.end_pos()) break;
        if ((targetStartOffset+static_cast<pos_t>(i)) <= scoreRange.begin_pos()) continue;

        const char targetBase(targetSeq[targetStartOffset+i]);

        if ((querySeq[i] != targetBase) || (querySeq[i] == 'N'))
        {
            if ((querySeq[i] == 'N') || (targetBase == 'N'))
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
    const std::string& targetSeq,
    const unsigned bestPos,
    SRAlignmentInfo& alignment)
{
    const unsigned querySize = querySeq.size();
    alignment.leftMismatches=0;
    alignment.homMismatches=0;
    alignment.rightMismatches=0;

    assert(bestPos+querySize <= targetSeq.size());

    for (unsigned i(0); i<querySize; i++)
    {
        if ((querySeq[i] != targetSeq[bestPos+i]) || (querySeq[i] == 'N'))
        {
            if (i<=alignment.leftSize)
            {
                alignment.leftMismatches += 1;
            }
            else if (i <= (alignment.leftSize+alignment.homSize))
            {
                alignment.homMismatches += 1;
            }
            else
            {
                alignment.rightMismatches += 1;
            }
        }
    }

    alignment.alignScore = querySize - (alignment.leftMismatches+ alignment.homMismatches+ alignment.rightMismatches);
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
    static const unsigned minFlankSize(16);
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
    const unsigned flankScoreSize,
    const std::string& querySeq,
    const qscore_snp& qualConvert,
    const uint8_t* queryQual,
    const std::string& targetSeq,
    const known_pos_range2& targetBpOffsetRange,
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
            << "querySeq:\n";
        printSeq(querySeq,oss);
        oss << '\n'
            << "targetSeq:\n";
        printSeq(targetSeq,oss);
        oss << '\n';
        BOOST_THROW_EXCEPTION(LogicException(oss.str()));
    }

    // set the scanning start & end to make sure the candidate windows overlapping the breakpoint
    const unsigned scanStart(std::max(0, static_cast<pos_t>(targetBpOffsetRange.begin_pos()) - static_cast<pos_t>(querySize) + 2));
    const unsigned scanEnd(std::max(0, std::min((static_cast<pos_t>(targetBpOffsetRange.end_pos())), static_cast<pos_t>(targetSize - querySize))));

    const known_pos_range2 scoreRange(targetBpOffsetRange.begin_pos()-static_cast<pos_t>(flankScoreSize),
                                      targetBpOffsetRange.end_pos()+static_cast<pos_t>(flankScoreSize));

#ifdef DEBUG_SRA
    static const std::string logtag("splitReadAligner: ");
    log_os << logtag << "query size = " << querySize << " target size = " << targetSize << '\n';
    log_os << logtag << "targetBeginPos = " << targetBpBeginPos << '\n';
    log_os << logtag << "scan start = " << scanStart << " scan end = " << scanEnd << '\n';
#endif
    if (scanEnd < scanStart)
    {
        std::ostringstream oss;
        oss << "ERROR: scanEnd < scanStart."
            << " scanEnd: " << scanEnd << " scanStart: " << scanStart
            << " querySize: " << querySize << " targetSize: " << targetSize << '\n'
            << "\ttargetRange: " << targetBpOffsetRange << '\n';
        BOOST_THROW_EXCEPTION(LogicException(oss.str()));
    }

    // do one high-speed pass to find the optimal alignment (in terms of lhood), then compute all the goodies later:
    float bestLnLhood(0);
    unsigned bestPos(0);
    {
        bool isBest(false);
        for (unsigned i = scanStart; i<= scanEnd; i++)
        {
            const float lnLhood(getLnLhood(querySeq, qualConvert, queryQual,
                                           targetSeq, i, scoreRange, isBest, bestLnLhood));

#ifdef DEBUG_SRA
            log_os << logtag << "scanning: " << i << " lhood: " << lnLhood << " bestLnLhood " << bestLnLhood << " isBest " << isBest << " bestPos " << bestPos << '\n';
#endif
            if ( (! isBest) || (lnLhood > bestLnLhood))
            {
                bestLnLhood = lnLhood;
                bestPos=i;
                isBest=true;
            }
        }
        assert(isBest);
    }

    // align query to itself to get the 'perfect' lhood:
    //alignment.perfectLnLhood = (getLnLhood(querySeq, qualConvert, queryQual,
    //                               querySeq.begin(), querySeq.end(), false, 0));

    assert(static_cast<pos_t>(bestPos) <= (targetBpOffsetRange.end_pos()+1));
    if (static_cast<pos_t>(bestPos) <= (targetBpOffsetRange.begin_pos()+1))
    {
        alignment.leftSize = static_cast<pos_t>(targetBpOffsetRange.begin_pos()+1) - bestPos;
    }
    else
    {
        alignment.leftSize = 0;
    }

    if (alignment.leftSize > querySize)
    {
        std::ostringstream oss;
        oss << "ERROR: Unexpected split read alignment outcome. "
            << " targetRange: " << targetBpOffsetRange << " bestPos: " << bestPos << " bestLnLhood: " << bestLnLhood << " querySize: " << querySize << " targetSize: " << targetSize << '\n'
            << "alignment: " << alignment << "\n"
            << "querySeq:\n";
        printSeq(querySeq,oss);
        oss << '\n'
            << "targetSeq:\n";
        printSeq(targetSeq,oss);
        oss << '\n';
        BOOST_THROW_EXCEPTION(LogicException(oss.str()));
    }
    alignment.homSize = std::min(querySize-alignment.leftSize,(static_cast<pos_t>(targetBpOffsetRange.end_pos()+1) - bestPos) - alignment.leftSize);

    if ((alignment.leftSize + alignment.homSize) < querySize)
    {
        alignment.rightSize = querySize - (alignment.leftSize + alignment.homSize);
    }
    else
    {
        alignment.rightSize = 0;
    }
    alignment.alignLnLhood = bestLnLhood;
    alignment.alignPos = bestPos;

    calculateAlignScore(querySeq, targetSeq, bestPos, alignment);

    // filtering the alignment and set evidence
    setEvidence(alignment);

#ifdef DEBUG_SRA
    log_os << logtag << "bestpos: " << bestPos << " final alignment\n" << alignment << "\n";
#endif
}
