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

#include <cassert>

#include <iostream>

#ifdef ALN_DEBUG
#include "blt_util/log.hh"
#include <iostream>
#endif


template <typename ScoreType>
std::ostream&
operator<<(std::ostream& os, JumpAlignmentResult<ScoreType>& alignment)
{
    os << "JumpAligner: score: " << alignment.score << "\n"
       << "\talign1: " << alignment.align1 << "\n"
       << "\talign2: " << alignment.align2 << "\n"
       << "\tjumpInsertSize " << alignment.jumpInsertSize << "\n"
       << "\tjumpRange " << alignment.jumpRange << "\n";
    return os;
}



// traceback:
template <typename ScoreType>
template <typename SymIter, typename MatrixType>
void
JumpAlignerBase<ScoreType>::
backTraceAlignment(
    const SymIter queryBegin, const SymIter queryEnd,
    const SymIter ref1Begin, const SymIter ref1End,
    const SymIter ref2Begin, const SymIter ref2End,
    const size_t querySize, const size_t ref1Size, const size_t ref2Size,
    const MatrixType& ptrMatrix1, const MatrixType& ptrMatrix2,
    const BackTrace<ScoreType>& btraceInput,
    JumpAlignmentResult<ScoreType>& result) const
{
    BackTrace<ScoreType> btrace(btraceInput);

    assert(btrace.isInit);
    assert(btrace.refBegin <= ref1Size+ref2Size);
    assert(btrace.queryBegin <= querySize);

    result.score = btrace.max;

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
        const MatrixType* ptrMatrixX(isRef1 ? &ptrMatrix1 : &ptrMatrix2 );
        const AlignState::index_t nextState(static_cast<AlignState::index_t>(ptrMatrixX->val(btrace.queryBegin,refXBegin).get(btrace.state)));

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
