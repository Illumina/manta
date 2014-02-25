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

//
// \author Chris Saunders
//

//#define ALN_DEBUG

#include "AlignerUtil.hh"

#include <cassert>

#include <iostream>

#ifdef ALN_DEBUG
#include "blt_util/log.hh"
#include <iostream>
#endif



template <typename ScoreType>
template <typename SymIter>
void
GlobalJumpIntronAligner<ScoreType>::
align(
    const SymIter queryBegin, const SymIter queryEnd,
    const SymIter ref1Begin, const SymIter ref1End,
    const SymIter ref2Begin, const SymIter ref2End,
    JumpAlignmentResult<ScoreType>& result) const
{
    result.clear();

    const AlignmentScores<ScoreType>& scores(this->getScores());

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
    // be soft-clipped and each base off the end will count as an 'offEdge' state:
    for (unsigned queryIndex(0); queryIndex<=querySize; queryIndex++)
    {
        ScoreVal& val((*thisSV)[queryIndex]);
        val.match = queryIndex * scores.offEdge;
        val.del = badVal;
        val.ins = badVal;
        val.jump = badVal;
        val.intron = badVal;
    }

    BackTrace<ScoreType> btrace;

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
                val.intron = badVal;
            }

            unsigned queryIndex(0);
            for (SymIter queryIter(queryBegin); queryIter != queryEnd; ++queryIter, ++queryIndex)
            {
                // update match
                ScoreVal& headScore((*thisSV)[queryIndex+1]);
                PtrVal& headPtr(_ptrMat1.val(queryIndex+1,ref1Index+1));
                {
                    const ScoreVal& sval((*prevSV)[queryIndex]);
                    headPtr.match = this->max3(
                                        headScore.match,
                                        sval.match,
                                        sval.del,
                                        sval.ins);
                    // Only can leave the intron (splice) state if the last two
                    // bases of the intron match the motif
                    if (*(ref1Iter-2)=='A' && *(ref1Iter-1)=='G')
                    {
                        if (headScore.match < sval.intron)
                        {
                            headScore.match = sval.intron;
                            headPtr.match = AlignState::SPLICE;
                        }
                    }

                    headScore.match += ((*queryIter==*ref1Iter) ? scores.match : scores.mismatch);
                }

                // update delete
                {
                    const ScoreVal& sval((*prevSV)[queryIndex+1]);
                    headPtr.del = this->max3(
                                      headScore.del,
                                      sval.match + scores.open,
                                      sval.del,
                                      sval.ins + scores.open);

                    headScore.del += scores.extend;
                    if (0==queryIndex) headScore.del += badVal;
                }

                // update insert
                {
                    const ScoreVal& sval((*thisSV)[queryIndex]);
                    headPtr.ins = this->max3(
                                      headScore.ins,
                                      sval.match + scores.open,
                                      sval.del + scores.open,
                                      sval.ins);

                    headScore.ins += scores.extend;
                    if (0==queryIndex) headScore.ins += badVal;
                }
                // update intron / splice
                {
                    const ScoreVal& sval((*prevSV)[queryIndex+1]);
                    headPtr.intron = AlignState::SPLICE;
                    headScore.intron = sval.intron;
                    // Only can enter the intron (splice) state if the first two
                    // bases of the intron match the motif
                    if (*(ref1Iter)=='G' && *(ref1Iter+1)=='T')
                    {
                        if (sval.match + _intronOpenScore > sval.intron)
                        {
                            headScore.intron = sval.match + _intronOpenScore;
                            headPtr.intron = AlignState::MATCH;
                        }
                    }
                }
                // update jump
                {
                    const ScoreVal& sval((*prevSV)[queryIndex+1]);
                    headPtr.jump = this->max4(
                                       headScore.jump,
                                       headScore.match + this->getJumpScore(),
                                       badVal,
                                       headScore.ins + this->getJumpScore(),
                                       sval.jump);
                }

#ifdef ALN_DEBUG
                log_os << "queryIdx refIdx ref1Idx: " << queryIndex+1 << " " << ref1Index+1 << " " << ref1Index+1 << "\n";
                log_os << headScore.match << ":" << headScore.del << ":" << headScore.ins << ":" << headScore.jump << ":" << headScore.intron << "/"
                       << static_cast<int>(headPtr.match) << static_cast<int>(headPtr.del) << static_cast<int>(headPtr.ins) << static_cast<int>(headPtr.jump) << "\n";
                log_os << "Q:" << *queryIter << " R:" << *ref1Iter << "\n";
#endif
            }
#ifdef ALN_DEBUG
            log_os << "\n";
#endif

            // get backtrace info:
            {
                const ScoreVal& sval((*thisSV)[querySize]);
                const ScoreType thisMax(sval.match);
                if (btrace.isInit && (thisMax<=btrace.max)) continue;
                btrace.max=thisMax;
                btrace.refBegin=ref1Index+1;
                btrace.queryBegin=querySize;
                btrace.isInit=true;
            }
        }
    }

    // in the backtrace start search, also allow for the case where the query falls-off the end of the reference:
    for (unsigned queryIndex(0); queryIndex<=querySize; queryIndex++)
    {
        const ScoreVal& sval((*thisSV)[queryIndex]);
        const ScoreType thisMax(sval.match + (querySize-queryIndex) * scores.offEdge);
        if (btrace.isInit && (thisMax<=btrace.max)) continue;
        btrace.max=thisMax;
        btrace.refBegin=ref1Size;
        btrace.queryBegin=queryIndex;
        btrace.isInit=true;
    }


    {
        for (unsigned queryIndex(0); queryIndex<=querySize; queryIndex++)
        {
            ScoreVal& val((*thisSV)[queryIndex]);
            val.match = queryIndex * scores.offEdge;
            val.del = badVal;
            val.ins = badVal;
            val.intron = badVal;
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
                val.intron = badVal;
            }

            unsigned queryIndex(0);
            for (SymIter queryIter(queryBegin); queryIter != queryEnd; ++queryIter, ++queryIndex)
            {
                // update match
                ScoreVal& headScore((*thisSV)[queryIndex+1]);
                PtrVal& headPtr(_ptrMat2.val(queryIndex+1,ref2Index+1));
                {
                    const ScoreVal& sval((*prevSV)[queryIndex]);
                    headPtr.match = this->max4(
                                        headScore.match,
                                        sval.match,
                                        sval.del,
                                        sval.ins,
                                        sval.jump);
                    // Only can leave the intron (splice) state if the last two
                    // bases of the intron match the motif
                    if (*(ref2Iter-2)=='A' && *(ref2Iter-1)=='G')
                    {
                        if (headScore.match < sval.intron)
                        {
                            headScore.match = sval.intron;
                            headPtr.match = AlignState::SPLICE;
                        }
                    }

                    headScore.match += ((*queryIter==*ref2Iter) ? scores.match : scores.mismatch);
                }

                // update delete
                {
                    const ScoreVal& sval((*prevSV)[queryIndex+1]);
                    headPtr.del = this->max3(
                                      headScore.del,
                                      sval.match + scores.open,
                                      sval.del,
                                      sval.ins + scores.open);

                    headScore.del += scores.extend;
                }

                // update insert
                {
                    const ScoreVal& sval((*thisSV)[queryIndex]);
                    headPtr.ins = this->max4(
                                      headScore.ins,
                                      sval.match + scores.open,
                                      sval.del + scores.open,
                                      sval.ins,
                                      sval.jump); // jump->ins moves get a pass on the gap-open penalty, to support mirco-insertions

                    headScore.ins += scores.extend;
                }
                // update intron
                {
                    const ScoreVal& sval((*prevSV)[queryIndex+1]);
                    headPtr.intron = AlignState::SPLICE;
                    headScore.intron = sval.intron;
                    // Only can enter the intron (splice) state if the first two
                    // bases of the intron match the motif
                    if (*(ref2Iter)=='G' && *(ref2Iter+1)=='T')
                    {
                        if (sval.match + _intronOpenScore > sval.intron)
                        {
                            headScore.intron = sval.match + _intronOpenScore;
                            headPtr.intron = AlignState::MATCH;
                        }
                    }
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
                if (btrace.isInit && (thisMax<=btrace.max)) continue;
                btrace.max=thisMax;
                btrace.refBegin=ref1Size+ref2Index+1;
                btrace.queryBegin=querySize;
                btrace.isInit=true;
            }
        }

    }

    // in the backtrace start search, also allow for the case where the query falls-off the end of the reference:
    for (unsigned queryIndex(0); queryIndex<=querySize; queryIndex++)
    {
        const ScoreVal& sval((*thisSV)[queryIndex]);
        const ScoreType thisMax(sval.match + (querySize-queryIndex) * scores.offEdge);
        if (btrace.isInit && (thisMax<=btrace.max)) continue;
        btrace.max=thisMax;
        btrace.refBegin=ref1Size+ref2Size;
        btrace.queryBegin=queryIndex;
        btrace.isInit=true;
    }

    assert(btrace.isInit);
    assert(btrace.refBegin <= ref1Size+ref2Size);
    assert(btrace.queryBegin <= querySize);

    result.score = btrace.max;

#ifdef ALN_DEBUG
    log_os << "bt-start queryIndex: " << btrace.queryBegin << " refIndex: " << btrace.refBegin << " state: " << AlignState::label(btrace.state) << " maxScore: " << btrace.max << "\n";
#endif

    // traceback:
    ALIGNPATH::path_t& apath1(result.align1.apath);
    ALIGNPATH::path_t& apath2(result.align2.apath);
    ALIGNPATH::path_segment ps;

    // add any trailing soft-clip if we go off the end of the reference:
    if (btrace.queryBegin < querySize)
    {
        ps.type = ALIGNPATH::SOFT_CLIP;
        ps.length = (querySize-btrace.queryBegin);
    }

    bool isRef2End(false);

    while ((btrace.queryBegin>0) && (btrace.refBegin>0))
    {
        if (isRef2End) break;
        const bool isRef1(btrace.refBegin<=ref1Size);
        ALIGNPATH::path_t& apath( isRef1 ? apath1 : apath2 );
        const unsigned refXBegin(btrace.refBegin - (isRef1 ? 0 : ref1Size));
        const PtrMat* ptrMatX(isRef1 ? &_ptrMat1 : &_ptrMat2 );
        const AlignState::index_t nextState(static_cast<AlignState::index_t>(ptrMatX->val(btrace.queryBegin,refXBegin).get(btrace.state)));

#ifdef ALN_DEBUG
        log_os << "bt-iter queryIndex: " << btrace.queryBegin
               << " refIndex: " << btrace.refBegin
               << " state: " << AlignState::label(btrace.state)
               << " next: " << AlignState::label(nextState)
               << "\n";
        log_os << "\tisref1: " << isRef1 << " refXBegin: " << refXBegin << "\n";
#endif

        if      (btrace.state==AlignState::MATCH)
        {
            if ((!isRef1) && (refXBegin==1) && (nextState==AlignState::MATCH)) isRef2End=true;

            AlignerUtil::updatePath(apath,ps,ALIGNPATH::MATCH);
            btrace.queryBegin--;
            btrace.refBegin--;
        }
        else if (btrace.state==AlignState::DELETE)
        {
            AlignerUtil::updatePath(apath,ps,ALIGNPATH::DELETE);
            btrace.refBegin--;
        }
        else if (btrace.state==AlignState::SPLICE)
        {
            AlignerUtil::updatePath(apath,ps,ALIGNPATH::SKIP);
            btrace.refBegin--;
        }
        else if (btrace.state==AlignState::INSERT)
        {
            AlignerUtil::updatePath(apath,ps,ALIGNPATH::INSERT);
            btrace.queryBegin--;
        }
        else if (btrace.state==AlignState::JUMP)
        {
            if (ps.type != ALIGNPATH::NONE)
            {
                assert(btrace.refBegin>=ref1Size);
                result.align2.beginPos = btrace.refBegin-ref1Size;
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
                if (nextState == AlignState::JUMP) btrace.refBegin--;
            }
        }
        else
        {
            assert(false && "Unknown align state");
        }
        btrace.state=nextState;
        ps.length++;
    }

    const bool isRef1(btrace.refBegin<ref1Size);
    ALIGNPATH::path_t& apath( isRef1 ? apath1 : apath2 );

    if (ps.type != ALIGNPATH::NONE) apath.push_back(ps);

    // soft-clip beginning of read if we fall off the end of the reference
    if (btrace.queryBegin!=0)
    {
        ps.type = ALIGNPATH::SOFT_CLIP;
        ps.length = btrace.queryBegin;
        apath.push_back(ps);
    }

    if (isRef1)
    {
        result.align1.beginPos = btrace.refBegin;
    }
    else
    {
        result.align2.beginPos = btrace.refBegin-ref1Size;
    }

    std::reverse(apath1.begin(),apath1.end());
    std::reverse(apath2.begin(),apath2.end());

    // figure out jumpRange:
    if (result.align1.isAligned() && result.align2.isAligned())
    {
        // find the distance over which ref1 and ref2 are equal following the start of the breakpoint
        SymIter ref1JumpIter(ref1Begin + result.align1.beginPos + apath_ref_length(apath1));
        SymIter ref2JumpIter(ref2Begin + result.align2.beginPos);
        while (true)
        {
            if (ref1JumpIter == ref1End) break;
            if (ref2JumpIter == ref2End) break;
            if ((*ref1JumpIter) != (*ref2JumpIter)) break;

            result.jumpRange++;
            ref1JumpIter++;
            ref2JumpIter++;
        }
    }

    // if true, output final cigars using seq match '=' and mismatch 'X' symbols:
    static const bool isOutputSeqMatch(true);

    if (isOutputSeqMatch)
    {
        apath_add_seqmatch(queryBegin, queryEnd, (ref1Begin+result.align1.beginPos), ref1End, apath1);

        const unsigned queryOffset = apath_read_length(apath1) + result.jumpInsertSize;
        apath_add_seqmatch(queryBegin + queryOffset, queryEnd, (ref2Begin+result.align2.beginPos), ref2End, apath2);
    }
}

