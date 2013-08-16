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
#include <iostream>
#endif


template <typename ScoreType>
std::ostream&
operator<<(std::ostream& os, JumpAlignmentResult<ScoreType>& alignment)
{
    os << "JumpAligner: score: " << alignment.score << "\n"
       << "\talign1: " << alignment.align1 << "\n"
       << "\talign2: " << alignment.align2 << "\n"
       << "\tjumpInsertSize " << alignment.jumpInsertSize << "\n";
    return os;
}



template <typename ScoreType>
template <typename SymIter>
void
GlobalJumpAligner<ScoreType>::
align(
    const SymIter queryBegin, const SymIter queryEnd,
    const SymIter ref1Begin, const SymIter ref1End,
    const SymIter ref2Begin, const SymIter ref2End,
    JumpAlignmentResult<ScoreType>& result) const
{
#ifdef ALN_DEBUG
    std::ostream& log_os(std::cerr);
#endif

    result.clear();

    const size_t querySize(std::distance(queryBegin, queryEnd));
    const size_t ref1Size(std::distance(ref1Begin, ref1End));
    const size_t ref2Size(std::distance(ref2Begin, ref2End));

    assert(0 != querySize);
    assert(0 != ref1Size);
    assert(0 != ref2Size);

    _score1.resize(querySize+1);
    _score2.resize(querySize+1);
    _ptrMat1.resize(querySize+1, ref1Size+1);
    _ptrMat2.resize(querySize+1, ref2Size+1);

    ScoreVec* thisSV(&_score1);
    ScoreVec* prevSV(&_score2);

    static const ScoreType badVal(-10000);

    // global alignment of query -- disallow start from insertion or deletion
    // state, query can 'fall-off' the end of a short reference, in which case it will
    // be soft-clipped and each base off the end will count as a mismatch:
    for (unsigned queryIndex(0); queryIndex<=querySize; queryIndex++)
    {
        ScoreVal& val((*thisSV)[queryIndex]);
        val.match = queryIndex * _scores.offEdge;
        val.del = badVal;
        val.ins = badVal;
        val.jump = badVal;
    }

    BackTrace<ScoreType> bt;

    {
        unsigned ref1Index(0);
        for (SymIter ref1Iter(ref1Begin); ref1Iter != ref1End; ++ref1Iter, ++ref1Index)
        {
            std::swap(thisSV,prevSV);

            {
                // disallow start from the insert or delete state:
                ScoreVal& val((*thisSV)[0]);
                val.match = 0;
                val.del = badVal;
                val.ins = badVal;
                val.jump = badVal;
            }

            unsigned queryIndex(0);
            for (SymIter queryIter(queryBegin); queryIter != queryEnd; ++queryIter, ++queryIndex)
            {
                // update match
                ScoreVal& headScore((*thisSV)[queryIndex+1]);
                PtrVal& headPtr(_ptrMat1.val(queryIndex+1,ref1Index+1));
                {
                    const ScoreVal& sval((*prevSV)[queryIndex]);
                    headPtr.match = max3(
                                        headScore.match,
                                        sval.match,
                                        sval.del,
                                        sval.ins);

                    headScore.match += ((*queryIter==*ref1Iter) ? _scores.match : _scores.mismatch);
                }

                // update delete
                {
                    const ScoreVal& sval((*prevSV)[queryIndex+1]);
                    headPtr.del = max3(
                                      headScore.del,
                                      sval.match + _scores.open,
                                      sval.del,
                                      sval.ins);

                    headScore.del += _scores.extend;
                    if (0==queryIndex) headScore.del += badVal;
                }

                // update insert
                {
                    const ScoreVal& sval((*thisSV)[queryIndex]);
                    headPtr.ins = max3(
                                      headScore.ins,
                                      sval.match + _scores.open,
                                      sval.del,
                                      sval.ins);

                    headScore.ins += _scores.extend;
                    if (0==queryIndex) headScore.ins += badVal;
                }

                // update jump
                {
                    const ScoreVal& sval((*prevSV)[queryIndex+1]);
                    headPtr.jump = max4(
                                       headScore.jump,
                                       headScore.match + _jumpScore,
                                       badVal,
                                       headScore.ins + _jumpScore,
                                       sval.jump);
                }

#ifdef ALN_DEBUG
                log_os << "queryIdx refIdx ref1Idx: " << queryIndex+1 << " " << ref1Index+1 << " " << ref1Index+1 << "\n";
                log_os << headScore.match << ":" << headScore.del << ":" << headScore.ins << ":" << headScore.jump << "/"
                       << static_cast<int>(headPtr.match) << static_cast<int>(headPtr.del) << static_cast<int>(headPtr.ins) << static_cast<int>(headPtr.jump) << "\n";
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
                bt.refBegin=ref1Index+1;
                bt.queryBegin=querySize;
                bt.isInit=true;
            }
        }
    }

    // in the backtrace start search, also allow for the case where the query falls-off the end of the reference:
    for (unsigned queryIndex(0); queryIndex<=querySize; queryIndex++)
    {
        const ScoreVal& sval((*thisSV)[queryIndex]);
        const ScoreType thisMax(sval.match + (querySize-queryIndex) * _scores.offEdge);
        if (bt.isInit && (thisMax<=bt.max)) continue;
        bt.max=thisMax;
        bt.refBegin=ref1Size;
        bt.queryBegin=queryIndex;
        bt.isInit=true;
    }


    {
        for (unsigned queryIndex(0); queryIndex<=querySize; queryIndex++)
        {
            ScoreVal& val((*thisSV)[queryIndex]);
            val.match = queryIndex * _scores.offEdge;
            val.del = badVal;
            val.ins = badVal;
            //val.jump = badVal; // preserve jump setting from last iteration of ref1
        }

        unsigned ref2Index(0);
        for (SymIter ref2Iter(ref2Begin); ref2Iter != ref2End; ++ref2Iter, ++ref2Index)
        {
            std::swap(thisSV,prevSV);

            {
                // disallow start from the insert or delete state:
                ScoreVal& val((*thisSV)[0]);
                val.match = 0;
                val.del = badVal;
                val.ins = badVal;
                val.jump = badVal;
            }

            unsigned queryIndex(0);
            for (SymIter queryIter(queryBegin); queryIter != queryEnd; ++queryIter, ++queryIndex)
            {
                // update match
                ScoreVal& headScore((*thisSV)[queryIndex+1]);
                PtrVal& headPtr(_ptrMat2.val(queryIndex+1,ref2Index+1));
                {
                    const ScoreVal& sval((*prevSV)[queryIndex]);
                    headPtr.match = max4(
                                        headScore.match,
                                        sval.match,
                                        sval.del,
                                        sval.ins,
                                        sval.jump);

                    headScore.match += ((*queryIter==*ref2Iter) ? _scores.match : _scores.mismatch);
                }

                // update delete
                {
                    const ScoreVal& sval((*prevSV)[queryIndex+1]);
                    headPtr.del = max3(
                                      headScore.del,
                                      sval.match + _scores.open,
                                      sval.del,
                                      sval.ins);

                    headScore.del += _scores.extend;
                }

                // update insert
                {
                    const ScoreVal& sval((*thisSV)[queryIndex]);
                    headPtr.ins = max4(
                                      headScore.ins,
                                      sval.match + _scores.open,
                                      sval.del,
                                      sval.ins,
                                      sval.jump + _scores.open);

                    headScore.ins += _scores.extend;
                }

                // update jump
                {
                    const ScoreVal& sval((*prevSV)[queryIndex+1]);
                    headPtr.jump = AlignState::JUMP;
                    headScore.jump = sval.jump;
                }

#ifdef ALN_DEBUG
                log_os << "queryIdx refIdx ref2Idx: " << queryIndex+1 << " " << ref1Size+ref2Index+1 << " " << ref2Index+1 << "\n";
                log_os << headScore.match << ":" << headScore.del << ":" << headScore.ins << ":" << headScore.jump << "/"
                       << static_cast<int>(headPtr.match) << static_cast<int>(headPtr.del) << static_cast<int>(headPtr.ins) << static_cast<int>(headPtr.jump) << "\n";
#endif
            }
#ifdef ALN_DEBUG
            log_os << "\n";
#endif

            // get backtrace start info:
            {
                const ScoreVal& sval((*thisSV)[querySize]);
                const ScoreType thisMax(sval.match);
                if (bt.isInit && (thisMax<=bt.max)) continue;
                bt.max=thisMax;
                bt.refBegin=ref1Size+ref2Index+1;
                bt.queryBegin=querySize;
                bt.isInit=true;
            }
        }

    }

    // in the backtrace start search, also allow for the case where the query falls-off the end of the reference:
    for (unsigned queryIndex(0); queryIndex<=querySize; queryIndex++)
    {
        const ScoreVal& sval((*thisSV)[queryIndex]);
        const ScoreType thisMax(sval.match + (querySize-queryIndex) * _scores.offEdge);
        if (bt.isInit && (thisMax<=bt.max)) continue;
        bt.max=thisMax;
        bt.refBegin=ref1Size+ref2Size;
        bt.queryBegin=queryIndex;
        bt.isInit=true;
    }

    assert(bt.isInit);
    assert(bt.refBegin <= ref1Size+ref2Size);
    assert(bt.queryBegin <= querySize);

    result.score = bt.max;

#ifdef ALN_DEBUG
    log_os << "bt-start queryIndex: " << bt.queryBegin << " refIndex: " << bt.refBegin << " state: " << AlignState::label(bt.state) << " maxScore: " << bt.max << "\n";
#endif

    // traceback:
    ALIGNPATH::path_t& apath1(result.align1.apath);
    ALIGNPATH::path_t& apath2(result.align2.apath);
    ALIGNPATH::path_segment ps;

    // add any trailing soft-clip if we go off the end of the reference:
    if (bt.queryBegin < querySize)
    {
        ps.type = ALIGNPATH::SOFT_CLIP;
        ps.length = (querySize-bt.queryBegin);
    }

    while ((bt.queryBegin>0) && (bt.refBegin>0))
    {
        const bool isRef1(bt.refBegin<=ref1Size);
        ALIGNPATH::path_t& apath( isRef1 ? apath1 : apath2 );
        const unsigned refXBegin(bt.refBegin - (isRef1 ? 0 : ref1Size));
        const PtrMat* ptrMatX(isRef1 ? &_ptrMat1 : &_ptrMat2 );
        const AlignState::index_t nextState(static_cast<AlignState::index_t>(ptrMatX->val(bt.queryBegin,refXBegin).get(bt.state)));

#ifdef ALN_DEBUG
        log_os << "bt-iter queryIndex: " << bt.queryBegin
               << " refIndex: " << bt.refBegin
               << " state: " << AlignState::label(bt.state)
               << " next: " << AlignState::label(nextState)
               << "\n";
        log_os << "\tisref1: " << isRef1 << " refXBegin: " << refXBegin << "\n";
#endif

        if      (bt.state==AlignState::MATCH)
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
        else if (bt.state==AlignState::JUMP)
        {
            if (ps.type != ALIGNPATH::NONE)
            {
                assert(bt.refBegin>=ref1Size);
                result.align2.beginPos = bt.refBegin-ref1Size;
                if (ps.type == ALIGNPATH::INSERT)
                {
                    result.jumpInsertSize += ps.length;
                    ps.type = ALIGNPATH::NONE;
                    ps.length = 0;
                }
                else
                {
                    AlignerUtil::updatePath(apath2,ps,ALIGNPATH::NONE);
                }
            }
            else
            {
                if (nextState == AlignState::JUMP) bt.refBegin--;
            }
        }
        else
        {
            assert(! "Unknown align state");
        }
        bt.state=nextState;
        ps.length++;
    }

    const bool isRef1(bt.refBegin<=ref1Size);
    ALIGNPATH::path_t& apath( isRef1 ? apath1 : apath2 );

    if (ps.type != ALIGNPATH::NONE) apath.push_back(ps);

    // soft-clip beginning of read if we fall off the end of the reference
    if (bt.queryBegin!=0)
    {
        ps.type = ALIGNPATH::SOFT_CLIP;
        ps.length = bt.queryBegin;
        apath.push_back(ps);
    }

    if (isRef1)
    {
        result.align1.beginPos = bt.refBegin;
    }
    else
    {
        result.align2.beginPos = bt.refBegin-ref1Size;
    }

    std::reverse(apath1.begin(),apath1.end());
    std::reverse(apath2.begin(),apath2.end());
}

