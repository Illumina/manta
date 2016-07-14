// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Manta - Structural Variant and Indel Caller
// Copyright (c) 2013-2016 Illumina, Inc.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//

///
/// \author Chris Saunders
///

#include <cassert>

#ifdef DEBUG_ALN_MATRIX
#include <boost/io/ios_state.hpp>

#include <iomanip>
#include <iostream>
#endif


#ifdef DEBUG_ALN_MATRIX
template <typename ScoreType>
template <typename SymIter, typename MatrixType, typename ScoreValType>
void
AlignerBase<ScoreType>::
dumpSingleRefTable(
    const SymIter refBegin, const SymIter refEnd,
    const size_t querySize,
    const MatrixType& ptrMatrix,
    const std::vector<std::vector<ScoreValType>>& storeScores,
    const char refSym,
    const AlignState::index_t sIndex,
    unsigned& storeIndex,
    std::ostream& os) const
{
    boost::io::ios_all_saver guard(os);

    auto printVal = [](
                    const ScoreType& val,
                    const char fromSym,
                    std::ostream& pos)
    {
        if (val<-900)
        {
            pos << " XX";
        }
        else
        {
            pos << std::setfill(' ') << std::setw(3) << val;
        }
        pos << fromSym;
    };

    auto printQueryRow = [&](
        const unsigned qrefIndex)
    {
        for (unsigned queryIndex(0); queryIndex <= querySize; ++queryIndex)
        {
            const auto& val(storeScores[storeIndex][queryIndex].getScore(sIndex));
            const char fromSym(AlignState::symbol(ptrMatrix.val(queryIndex, qrefIndex).getStatePtr(sIndex)));
            printVal(val, fromSym, os);
        }
        os << "\n";
    };

    os << "# - ";
    printQueryRow(0);
    unsigned refIndex(0);
    for (SymIter refIter(refBegin); refIter != refEnd; ++refIter, ++refIndex)
    {
        os << refSym << " " << *refIter << " ";
        storeIndex++;
        printQueryRow(refIndex+1);
    }
}
#endif



template <typename ScoreType>
ScoreType
AlignerBase<ScoreType>::
getPathScore(
    const ALIGNPATH::path_t& apath,
    const bool isScoreOffEdge) const
{
    using namespace ALIGNPATH;

    ScoreType val(0);

    // note that the intent of this function was to replicate the underling aligner and thus
    // not penalize insert-delete transitions, however as written the open score is added twice for
    // such an event. This turns out to perform much better for the performance of this tool as
    // a variant 'arm' validator. so the unintended code is staying in place for now.
    //
    /// TODO: reevaluate policy for insertion-deletion state transition score

    for (const path_segment& ps : apath)
    {
        bool isIndel(false); // placement of isIndel inside of this loop is the 'bug'
        switch (ps.type)
        {
        case MATCH:
            assert(false && "Unexpected MATCH segment"); // if MATCH segments exist, then you're using the wrong type of CIGAR for this function
            break;
        case SEQ_MATCH:
            val += (_scores.match * ps.length);
            isIndel = false;
            break;
        case SEQ_MISMATCH:
            val += (_scores.mismatch * ps.length);
            isIndel = false;
            break;
        case INSERT:
        case DELETE:
            if (! isIndel) val += _scores.open;
            val += (_scores.extend * ps.length);
            isIndel = true;
            break;
        case SOFT_CLIP:
            if (isScoreOffEdge) val += (_scores.offEdge * ps.length);
            isIndel = false;
            break;
        default:
            break;
        }
    }
    return val;
}



template <typename ScoreType>
ScoreType
AlignerBase<ScoreType>::
getMaxPathScore(
    const ALIGNPATH::path_t& apath,
    unsigned& maxReadOffset,
    unsigned& maxRefOffset,
    const bool isScoreOffEdge) const
{
    using namespace ALIGNPATH;

    ScoreType val(0);
    unsigned readOffset(0);
    unsigned refOffset(0);

    ScoreType maxVal(0);
    maxReadOffset=0;
    maxRefOffset=0;

    for (const path_segment& ps : apath)
    {
        bool isIndel(false); // unintended 'bug' with positive results, see TODO note above
        switch (ps.type)
        {
        case MATCH:
            assert(false && "Unexpected MATCH segment"); // if MATCH segments exist, then you're using the wrong type of CIGAR for this function
            break;
        case SEQ_MATCH:
            val += (_scores.match * ps.length);
            readOffset += ps.length;
            refOffset += ps.length;
            isIndel = false;
            break;
        case SEQ_MISMATCH:
            val += (_scores.mismatch * ps.length);
            readOffset += ps.length;
            refOffset += ps.length;
            isIndel = false;
            break;
        case INSERT:
            if (! isIndel) val += _scores.open;
            val += (_scores.extend * ps.length);
            readOffset += ps.length;
            isIndel = true;
            break;
        case DELETE:
            if (! isIndel) val += _scores.open;
            val += (_scores.extend * ps.length);
            refOffset += ps.length;
            isIndel = true;
            break;
        case SOFT_CLIP:
            if (isScoreOffEdge) val += (_scores.offEdge * ps.length);
            readOffset += ps.length;
            isIndel = false;
            break;
        default:
            break;
        }

        if (val>maxVal)
        {
            maxVal = val;
            maxReadOffset = readOffset;
            maxRefOffset = refOffset;
        }
    }
    return maxVal;
}
