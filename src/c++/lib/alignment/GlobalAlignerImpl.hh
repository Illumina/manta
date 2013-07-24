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

#include <limits>

#define ALN_DEBUG


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
    _scoreMat.resize(size1+1, size2+1);
    _ptrMat.resize(size1+1, size2+1);

    static const ScoreType badVal(-10000);

    // global alignment of seq1 -- disallow partial alignment:
    for (unsigned index1(0); index1<=size1; index1++)
    {
        ScoreVal& val(_scoreMat.val(index1,0));
        val.match = badVal;
        val.del = badVal;
        val.ins = badVal;
    }

    // disallow start from the delete state:
    for (unsigned index2(0); index2<=size2; index2++)
    {
        ScoreVal& val(_scoreMat.val(0,index2));
        val.match = 0;
        val.del = badVal;
        val.ins = 0;
    }

    unsigned index1(0);
    for (SymIter iter1(begin1); iter1 != end1; ++iter1, ++index1)
    {
        unsigned index2(0);
        for (SymIter iter2(begin2); iter2 != end2; ++iter2, ++index2)
        {
#ifdef ALN_DEBUG
            log_os << "i1i2: " << index1 << " " << index2 << "\n";
#endif
            // update match matrix
            ScoreVal& headScore(_scoreMat.val(index1+1,index2+1));
            PtrVal& headPtr(_ptrMat.val(index1+1,index2+1));
            {
                const ScoreVal& sval(_scoreMat.val(index1,index2));
                headPtr.match = max3(
                                    headScore.match,
                                    sval.match,
                                    sval.del + _scores.open,
                                    sval.ins + _scores.open);

                headScore.match += ((*iter1==*iter2) ? _scores.match : _scores.mismatch);
            }

            // update delete
            {
                const ScoreVal& sval(_scoreMat.val(index1+1,index2));
                headPtr.del = max3(
                                  headScore.del,
                                  sval.match,
                                  sval.del,
                                  sval.ins + _scores.open);

                headScore.del += _scores.extend;
            }

            // update insert
            {
                const ScoreVal& sval(_scoreMat.val(index1,index2+1));
                headPtr.ins = max3(
                                  headScore.ins,
                                  sval.match,
                                  sval.del + _scores.open,
                                  sval.ins);

                headScore.ins += _scores.extend;
            }

#ifdef ALN_DEBUG
            log_os << headScore.match << ":" << headScore.del << ":" << headScore.ins << "/"
                   << static_cast<int>(headPtr.match) << static_cast<int>(headPtr.del) << static_cast<int>(headPtr.ins) << "\n";
#endif
        }
#ifdef ALN_DEBUG
        log_os << "\n";
#endif
    }

    ScoreType max(0);
    matrix_t whichMatrix(MATCH);
    unsigned start1(size1+1);
    unsigned start2(0);

    for (unsigned index2(0); index2<=size2; index2++)
    {
        const ScoreVal& val(_scoreMat.val(size1+1, index2));
        ScoreType thisMax(val.match);
        matrix_t thisMatrix(MATCH);
        if (val.match <= val.ins)
        {
            thisMax=val.ins;
            thisMatrix=INSERT;
        }

        if ((index2==0) || (thisMax>max))
        {
            max=thisMax;
            start2=index2;
            whichMatrix=thisMatrix;
        }
    }

    result.score = max;

#ifdef ALN_DEBUG
    log_os << "start cell = " << start1 << " " << start2 << " " << whichMatrix << " " << max << "\n";
#endif

    // traceback:
    ALIGNPATH::path_segment ps;
    while ((start1>0) && (start2>0))
    {
        const matrix_t nextMatrix(static_cast<matrix_t>(_ptrMat.val(start1,start2).get(whichMatrix)));

        if (whichMatrix==MATCH)
        {
            updatePath(result.apath,ps,ALIGNPATH::MATCH);
            start1--;
            start2--;
        }
        else if (whichMatrix==DELETE)
        {
            updatePath(result.apath,ps,ALIGNPATH::DELETE);
            start2--;
        }
        else
        {
            updatePath(result.apath,ps,ALIGNPATH::INSERT);
            start1--;
        }
        ps.length++;
        whichMatrix=nextMatrix;
    }

    if (ps.type != ALIGNPATH::NONE) result.apath.push_back(ps);

    result.alignStart = start2;
    std::reverse(result.apath.begin(),result.apath.end());
}

