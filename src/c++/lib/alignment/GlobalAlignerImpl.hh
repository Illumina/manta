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

#ifdef ALN_DEBUG
#include <iostream>
#endif



template <typename ScoreType>
template <typename SymIter>
void
GlobalAligner<ScoreType>::
align(
    const SymIter queryBegin, const SymIter queryEnd,
    const SymIter refBegin, const SymIter refEnd,
    AlignmentResult<ScoreType>& result)
{
#ifdef ALN_DEBUG
    std::ostream& log_os(std::cerr);
#endif

    result.clear();

    const size_t querySize(std::distance(queryBegin, queryEnd));
    const size_t refSize(std::distance(refBegin, refEnd));

    assert(0 != querySize);
    assert(0 != refSize);

    _scoreMat.resize(querySize+1, refSize+1);
    _ptrMat.resize(querySize+1, refSize+1);

    static const ScoreType badVal(-10000);

    // global alignment of seq1 -- disallow start from insertion or deletion
    // state, seq1 can 'fall-off' the end of a short reference, in which case it will
    // be soft-clipped and each base off the end will count as a mismatch:
    for (unsigned queryIndex(0); queryIndex<=querySize; queryIndex++)
    {
        ScoreVal& val(_scoreMat.val(queryIndex,0));
        val.match = queryIndex * _scores.mismatch;
        val.del = badVal;
        val.ins = badVal;
    }

    // disallow start from the insert or delete state:
    for (unsigned refIndex(0); refIndex<=refSize; refIndex++)
    {
        ScoreVal& val(_scoreMat.val(0,refIndex));
        val.match = 0;
        val.del = badVal;
        val.ins = badVal;
    }

    {
        unsigned queryIndex(0);
        for (SymIter queryIter(queryBegin); queryIter != queryEnd; ++queryIter, ++queryIndex)
        {
            unsigned refIndex(0);
            for (SymIter refIter(refBegin); refIter != refEnd; ++refIter, ++refIndex)
            {
                // update match matrix
                ScoreVal& headScore(_scoreMat.val(queryIndex+1,refIndex+1));
                PtrVal& headPtr(_ptrMat.val(queryIndex+1,refIndex+1));
                {
                    const ScoreVal& sval(_scoreMat.val(queryIndex,refIndex));
                    headPtr.match = max3(
                                        headScore.match,
                                        sval.match,
                                        sval.del,
                                        sval.ins);

                    headScore.match += ((*queryIter==*refIter) ? _scores.match : _scores.mismatch);
                }

                // update delete
                {
                    const ScoreVal& sval(_scoreMat.val(queryIndex+1,refIndex));
                    headPtr.del = max3(
                                      headScore.del,
                                      sval.match + _scores.open,
                                      sval.del,
                                      sval.ins);

                    headScore.del += _scores.extend;
                }

                // update insert
                {
                    const ScoreVal& sval(_scoreMat.val(queryIndex,refIndex+1));
                    headPtr.ins = max3(
                                      headScore.ins,
                                      sval.match + _scores.open,
                                      sval.del,
                                      sval.ins);

                    headScore.ins += _scores.extend;
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
        }
    }

    //
    // find starting point for backtrace:
    //
    ScoreType max(0);
    AlignState::index_t whichMatrix(AlignState::MATCH);
    unsigned queryStart(0),refStart(0);
    bool isInit(false);

    for (unsigned refIndex(0); refIndex<=refSize; refIndex++)
    {
        const ScoreVal& val(_scoreMat.val(querySize, refIndex));
        const ScoreType thisMax(val.match);
        if(isInit && (thisMax<=max)) continue;
        max=thisMax;
        refStart=refIndex;
        queryStart=querySize;
        isInit=true;
    }

    // also allow for the case where seq1 falls-off the end of the reference:
    for (unsigned queryIndex(0); queryIndex<=querySize; queryIndex++)
    {
        const ScoreVal& val(_scoreMat.val(queryIndex, refSize));
        const ScoreType thisMax(val.match + (querySize-queryIndex) * _scores.mismatch);
        if(isInit && (thisMax<=max)) continue;
        max=thisMax;
        refStart=refSize;
        queryStart=queryIndex;
        isInit=true;
    }

    assert(isInit);
    assert(refStart <= refSize);
    assert(queryStart <= querySize);

    result.score = max;

#ifdef ALN_DEBUG
    log_os << "bt-start queryIndex: " << queryStart << " refIndex: " << refStart << " state: " << AlignState::label(whichMatrix) << " maxScore: " << max << "\n";
#endif

    // traceback:
    ALIGNPATH::path_t& apath(result.align.apath);
    ALIGNPATH::path_segment ps;

    // add any trailing soft-clip if we go off the end of the reference:
    if(queryStart < querySize)
    {
        ps.type = ALIGNPATH::SOFT_CLIP;
        ps.length = (querySize-queryStart);
    }

    while ((queryStart>0) && (refStart>0))
    {
        const AlignState::index_t nextMatrix(static_cast<AlignState::index_t>(_ptrMat.val(queryStart,refStart).get(whichMatrix)));

        if (whichMatrix==AlignState::MATCH)
        {
            AlignerUtil::updatePath(apath,ps,ALIGNPATH::MATCH);
            queryStart--;
            refStart--;
        }
        else if (whichMatrix==AlignState::DELETE)
        {
            AlignerUtil::updatePath(apath,ps,ALIGNPATH::DELETE);
            refStart--;
        }
        else if (whichMatrix==AlignState::INSERT)
        {
            AlignerUtil::updatePath(apath,ps,ALIGNPATH::INSERT);
            queryStart--;
        }
        else
        {
            assert(! "Unknown align state");
        }
        whichMatrix=nextMatrix;
        ps.length++;
    }

    if (ps.type != ALIGNPATH::NONE) apath.push_back(ps);

    // soft-clip beginning of read if we fall off the end of the reference
    if(queryStart!=0)
    {
        ps.type = ALIGNPATH::SOFT_CLIP;
        ps.length = queryStart;
        apath.push_back(ps);
    }

    result.align.alignStart = refStart;
    std::reverse(apath.begin(),apath.end());
}

