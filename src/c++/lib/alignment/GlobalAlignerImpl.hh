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


#ifdef ALN_DEBUG
#include <iostream>
#endif



template <typename ScoreType>
template <typename SymIter>
void
GlobalAligner<ScoreType>::
align(
    const SymIter begin1, const SymIter end1,
    const SymIter begin2, const SymIter end2,
    AlignmentResult<ScoreType>& result)
{
#ifdef ALN_DEBUG
    std::ostream& log_os(std::cerr);
#endif

    result.clear();

    const size_t size1(std::distance(begin1, end1));
    const size_t size2(std::distance(begin2, end2));

    assert(0 != size1);
    assert(0 != size2);

    _scoreMat.resize(size1+1, size2+1);
    _ptrMat.resize(size1+1, size2+1);

    static const ScoreType badVal(-10000);

    // global alignment of seq1 -- disallow start from insertion or deletion
    // state, seq1 can 'fall-off' the end of a short reference, in which case it will
    // be soft-clipped and each base off the end will count as a mismatch:
    for (unsigned index1(0); index1<=size1; index1++)
    {
        ScoreVal& val(_scoreMat.val(index1,0));
        val.match = index1 * _scores.mismatch;
        val.del = badVal;
        val.ins = badVal;
    }

    // disallow start from the insert or delete state:
    for (unsigned index2(0); index2<=size2; index2++)
    {
        ScoreVal& val(_scoreMat.val(0,index2));
        val.match = 0;
        val.del = badVal;
        val.ins = badVal;
    }

    {
        unsigned index1(0);
        for (SymIter iter1(begin1); iter1 != end1; ++iter1, ++index1)
        {
            unsigned index2(0);
            for (SymIter iter2(begin2); iter2 != end2; ++iter2, ++index2)
            {
                // update match matrix
                ScoreVal& headScore(_scoreMat.val(index1+1,index2+1));
                PtrVal& headPtr(_ptrMat.val(index1+1,index2+1));
                {
                    const ScoreVal& sval(_scoreMat.val(index1,index2));
                    headPtr.match = max3(
                                        headScore.match,
                                        sval.match,
                                        sval.del,
                                        sval.ins);

                    headScore.match += ((*iter1==*iter2) ? _scores.match : _scores.mismatch);
                }

                // update delete
                {
                    const ScoreVal& sval(_scoreMat.val(index1+1,index2));
                    headPtr.del = max3(
                                      headScore.del,
                                      sval.match + _scores.open,
                                      sval.del,
                                      sval.ins);

                    headScore.del += _scores.extend;
                }

                // update insert
                {
                    const ScoreVal& sval(_scoreMat.val(index1,index2+1));
                    headPtr.ins = max3(
                                      headScore.ins,
                                      sval.match + _scores.open,
                                      sval.del,
                                      sval.ins);

                    headScore.ins += _scores.extend;
                }

#ifdef ALN_DEBUG
                log_os << "i1i2: " << index1+1 << " " << index2+1 << "\n";
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
    unsigned start1(0),start2(0);
    bool isInit(false);

    for (unsigned index2(0); index2<=size2; index2++)
    {
        const ScoreVal& val(_scoreMat.val(size1, index2));
        const ScoreType thisMax(val.match);
        if(isInit && (thisMax<=max)) continue;
        max=thisMax;
        start2=index2;
        start1=size1;
        isInit=true;
    }

    // also allow for the case where seq1 falls-off the end of the reference:
    for (unsigned index1(0); index1<=size1; index1++)
    {
        const ScoreVal& val(_scoreMat.val(index1, size2));
        const ScoreType thisMax(val.match + (size1-index1) * _scores.mismatch);
        if(isInit && (thisMax<=max)) continue;
        max=thisMax;
        start2=size2;
        start1=index1;
        isInit=true;
    }

    assert(isInit);
    assert(start2 <= size2);
    assert(start1 <= size1);

    result.score = max;

#ifdef ALN_DEBUG
    log_os << "bt-start queryIndex: " << start1 << " refIndex: " << start2 << " state: " << AlignState::label(whichMatrix) << " maxScore: " << max << "\n";
#endif

    // traceback:
    ALIGNPATH::path_segment ps;

    // add any trailing soft-clip if we go off the end of the reference:
    if(start1 < size1)
    {
        ps.type = ALIGNPATH::SOFT_CLIP;
        ps.length = (size1-start1);
    }

    while ((start1>0) && (start2>0))
    {
        const AlignState::index_t nextMatrix(static_cast<AlignState::index_t>(_ptrMat.val(start1,start2).get(whichMatrix)));

        if (whichMatrix==AlignState::MATCH)
        {
            updatePath(result.apath,ps,ALIGNPATH::MATCH);
            start1--;
            start2--;
        }
        else if (whichMatrix==AlignState::DELETE)
        {
            updatePath(result.apath,ps,ALIGNPATH::DELETE);
            start2--;
        }
        else if (whichMatrix==AlignState::INSERT)
        {
            updatePath(result.apath,ps,ALIGNPATH::INSERT);
            start1--;
        }
        else
        {
            assert(! "Unknown align state");
        }
        whichMatrix=nextMatrix;
        ps.length++;
    }

    if (ps.type != ALIGNPATH::NONE) result.apath.push_back(ps);

    // soft-clip beginning of read if we fall off the end of the reference
    if(start1!=0)
    {
        ps.type = ALIGNPATH::SOFT_CLIP;
        ps.length = start1;
        result.apath.push_back(ps);
    }

    result.alignStart = start2;
    std::reverse(result.apath.begin(),result.apath.end());
}

