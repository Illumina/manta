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
// <https://github.com/downloads/sequencing/licenses/>.
//

/// derived from ELAND implementation by Tony Cox

//#define ALN_DEBUG

#include "AlignerUtil.hh"

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

    BackTrace<ScoreType> bt;

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
                if (bt.isInit && (thisMax<=bt.max)) continue;
                bt.max=thisMax;
                bt.refBegin=refIndex+1;
                bt.queryBegin=querySize;
                bt.isInit=true;
            }
        }
    }

    // also allow for the case where query falls-off the end of the reference:
    for (unsigned queryIndex(0); queryIndex<=querySize; queryIndex++)
    {
        const ScoreVal& sval((*thisSV)[queryIndex]);
        const ScoreType thisMax(sval.match + (querySize-queryIndex) * scores.offEdge);
        if (bt.isInit && (thisMax<=bt.max)) continue;
        bt.max=thisMax;
        bt.refBegin=refSize;
        bt.queryBegin=queryIndex;
        bt.isInit=true;
    }

    assert(bt.isInit);
    assert(bt.refBegin <= refSize);
    assert(bt.queryBegin <= querySize);

    result.score = bt.max;

#ifdef ALN_DEBUG
    log_os << "bt-start queryIndex: " << queryBegin << " refIndex: " << refBegin << " state: " << AlignState::label(bt.state) << " maxScore: " << bt.max << "\n";
#endif

    // traceback:
    ALIGNPATH::path_t& apath(result.align.apath);
    ALIGNPATH::path_segment ps;

    // add any trailing soft-clip if we go off the end of the reference:
    if (bt.queryBegin < querySize)
    {
        ps.type = ALIGNPATH::SOFT_CLIP;
        ps.length = (querySize-bt.queryBegin);
    }

    while ((bt.queryBegin>0) && (bt.refBegin>0))
    {
        const AlignState::index_t nextMatrix(static_cast<AlignState::index_t>(_ptrMat.val(bt.queryBegin,bt.refBegin).get(bt.state)));

        if (bt.state==AlignState::MATCH)
        {
            AlignerUtil::updatePath(apath,ps,ALIGNPATH::MATCH);
            bt.queryBegin--;
            bt.refBegin--;
        }
        else if (bt.state==AlignState::DELETE)
        {
            AlignerUtil::updatePath(apath,ps,ALIGNPATH::DELETE);
            bt.refBegin--;
        }
        else if (bt.state==AlignState::INSERT)
        {
            AlignerUtil::updatePath(apath,ps,ALIGNPATH::INSERT);
            bt.queryBegin--;
        }
        else
        {
            assert(! "Unknown align state");
        }
        bt.state=nextMatrix;
        ps.length++;
    }

    if (ps.type != ALIGNPATH::NONE) apath.push_back(ps);

    // soft-clip beginning of read if we fall off the end of the reference
    if (bt.queryBegin!=0)
    {
        ps.type = ALIGNPATH::SOFT_CLIP;
        ps.length = bt.queryBegin;
        apath.push_back(ps);
    }

    result.align.beginPos = bt.refBegin;
    std::reverse(apath.begin(),apath.end());

    // if true, output final cigars using seq match '=' and mismatch 'X' symbols:
    static const bool isOutputSeqMatch(true);

    if (isOutputSeqMatch)
    {
        apath_add_seqmatch(queryBegin, queryEnd, (refBegin+result.align.beginPos), refEnd, apath);
    }
}

