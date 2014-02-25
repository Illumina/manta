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

/// derived from ELAND implementation by Tony Cox

//#define ALN_DEBUG

#include <cassert>

#include <iostream>

#ifdef ALN_DEBUG
#include "blt_util/log.hh"
#include <iostream>
#endif



template <typename ScoreType>
std::ostream&
operator<<(std::ostream& os, AlignmentResult<ScoreType>& alignment)
{
    os << "AlignerResult: score: " << alignment.score << "\n"
       << "\talignment: " << alignment.align << "\n";
    return os;
}



// traceback:
template <typename ScoreType>
template <typename SymIter>
void
GlobalAligner<ScoreType>::
backTraceAlignment(
    const SymIter queryBegin, const SymIter queryEnd,
    const SymIter refBegin, const SymIter refEnd,
    const size_t querySize,
    const BackTrace<ScoreType>& btraceInput,
    Alignment& alignment) const
{
    BackTrace<ScoreType> btrace(btraceInput);

    // traceback:
    ALIGNPATH::path_t& apath(alignment.apath);
    ALIGNPATH::path_segment ps;

    // add any trailing soft-clip if we go off the end of the reference:
    if (btrace.queryBegin < querySize)
    {
        ps.type = ALIGNPATH::SOFT_CLIP;
        ps.length = (querySize-btrace.queryBegin);
    }

    while ((btrace.queryBegin>0) && (btrace.refBegin>0))
    {
        const AlignState::index_t nextMatrix(static_cast<AlignState::index_t>(_ptrMat.val(btrace.queryBegin,btrace.refBegin).get(btrace.state)));

        if (btrace.state==AlignState::MATCH)
        {
            AlignerUtil::updatePath(apath,ps,ALIGNPATH::MATCH);
            btrace.queryBegin--;
            btrace.refBegin--;
        }
        else if (btrace.state==AlignState::DELETE)
        {
            AlignerUtil::updatePath(apath,ps,ALIGNPATH::DELETE);
            btrace.refBegin--;
        }
        else if (btrace.state==AlignState::INSERT)
        {
            AlignerUtil::updatePath(apath,ps,ALIGNPATH::INSERT);
            btrace.queryBegin--;
        }
        else
        {
            assert(false && "Unknown align state");
        }
        btrace.state=nextMatrix;
        ps.length++;
    }

    if (ps.type != ALIGNPATH::NONE) apath.push_back(ps);

    // soft-clip beginning of read if we fall off the end of the reference
    if (btrace.queryBegin!=0)
    {
        ps.type = ALIGNPATH::SOFT_CLIP;
        ps.length = btrace.queryBegin;
        apath.push_back(ps);
    }

    alignment.beginPos = btrace.refBegin;
    std::reverse(apath.begin(),apath.end());

    // if true, output final cigars using seq match '=' and mismatch 'X' symbols:
    static const bool isOutputSeqMatch(true);

    if (isOutputSeqMatch)
    {
        apath_add_seqmatch(queryBegin, queryEnd, (refBegin+alignment.beginPos), refEnd, apath);
    }
}



template <typename ScoreType>
template <typename SymIter>
void
GlobalAligner<ScoreType>::
align(
    const SymIter queryBegin, const SymIter queryEnd,
    const SymIter refBegin, const SymIter refEnd,
    AlignmentResult<ScoreType>& result) const
{
    result.clear();

    const AlignmentScores<ScoreType>& scores(this->getScores());

    const size_t querySize(std::distance(queryBegin, queryEnd));
    const size_t refSize(std::distance(refBegin, refEnd));

    assert(0 != querySize);
    assert(0 != refSize);

    _score1.resize(querySize+1);
    _score2.resize(querySize+1);
    _ptrMat.resize(querySize+1, refSize+1);

    ScoreVec* thisSV(&_score1);
    ScoreVec* prevSV(&_score2);

    static const ScoreType badVal(-10000);

    // global alignment of query -- disallow start from insertion or deletion
    // state, query can 'fall-off' the end of a short reference, in which case it will
    // be soft-clipped and each base off the end will be scored as offEdge:
    for (unsigned queryIndex(0); queryIndex<=querySize; queryIndex++)
    {
        ScoreVal& val((*thisSV)[queryIndex]);
        val.match = queryIndex * scores.offEdge;
        val.del = badVal;
        val.ins = badVal;
    }

    BackTrace<ScoreType> btrace;

    {
        unsigned refIndex(0);
        for (SymIter refIter(refBegin); refIter != refEnd; ++refIter, ++refIndex)
        {
            std::swap(thisSV,prevSV);

            {
                // disallow start from the insert or delete state:
                ScoreVal& val((*thisSV)[0]);
                val.match = 0;
                val.del = badVal;
                val.ins = badVal;
            }

            unsigned queryIndex(0);
            for (SymIter queryIter(queryBegin); queryIter != queryEnd; ++queryIter, ++queryIndex)
            {
                // update match
                ScoreVal& headScore((*thisSV)[queryIndex+1]);
                PtrVal& headPtr(_ptrMat.val(queryIndex+1,refIndex+1));
                {
                    const ScoreVal& sval((*prevSV)[queryIndex]);
                    headPtr.match = this->max3(
                                        headScore.match,
                                        sval.match,
                                        sval.del,
                                        sval.ins);

                    headScore.match += ((*queryIter==*refIter) ? scores.match : scores.mismatch);
                }

                // update delete
                {
                    const ScoreVal& sval((*prevSV)[queryIndex+1]);
                    headPtr.del = this->max3(
                                      headScore.del,
                                      sval.match + scores.open,
                                      sval.del,
                                      sval.ins);

                    headScore.del += scores.extend;
                    if (0==queryIndex) headScore.del += badVal;
                }

                // update insert
                {
                    const ScoreVal& sval((*thisSV)[queryIndex]);
                    headPtr.ins = this->max3(
                                      headScore.ins,
                                      sval.match + scores.open,
                                      sval.del,
                                      sval.ins);

                    headScore.ins += scores.extend;
                    if (0==queryIndex) headScore.ins += badVal;
                }

#ifdef ALN_DEBUG
                log_os << "i1i2: " << queryIndex+1 << " " << refIndex+1 << "\n";
                log_os << headScore.match << ":" << headScore.del << ":" << headScore.ins << "/"
                       << static_cast<int>(headPtr.match) << static_cast<int>(headPtr.del) << static_cast<int>(headPtr.ins) << "\n";
#endif
            }
#ifdef ALN_DEBUG
            log_os << "\n";
#endif

            // get backtrace info:
            {
                const ScoreVal& sval((*thisSV)[querySize]);
                const ScoreType thisMax(sval.match);
                updateBacktrace(thisMax,refIndex+1,querySize,btrace);
            }
        }
    }

    // also allow for the case where query falls-off the end of the reference:
    for (unsigned queryIndex(0); queryIndex<=querySize; queryIndex++)
    {
        const ScoreVal& sval((*thisSV)[queryIndex]);
        const ScoreType thisMax(sval.match + (querySize-queryIndex) * scores.offEdge);
        updateBacktrace(thisMax,refSize,queryIndex,btrace);
    }

    assert(btrace.isInit);
    assert(btrace.refBegin <= refSize);
    assert(btrace.queryBegin <= querySize);


#ifdef ALN_DEBUG
    log_os << "btrace-start queryIndex: " << queryBegin << " refIndex: " << refBegin << " state: " << AlignState::label(btrace.state) << " maxScore: " << btrace.max << "\n";
#endif

    result.score = btrace.max;
    backTraceAlignment(queryBegin, queryEnd, refBegin, refEnd,
        querySize, btrace, result.align);
}

