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

#pragma once

#include "AlignerBase.hh"
#include "Alignment.hh"

#include "blt_util/basic_matrix.hh"

#include <iosfwd>


template <typename ScoreType>
struct AlignmentResult
{
    AlignmentResult()
    {
        clear();
    }

    void
    clear()
    {
        score=0;
        align.clear();
    }

    ScoreType score;
    Alignment align;
};


template <typename ScoreType>
std::ostream&
operator<<(std::ostream& os, AlignmentResult<ScoreType>& alignment);



/// \brief Implementation of global alignment with affine gap costs
///
/// alignment outputs start positions and CIGAR-style alignment
/// expression of query to reference. Alignment of
/// query is global -- query can "fall-off" either end of the reference,
/// in this case, each unaligned position is given an "off-edge" score and
/// the base is soft-clipped in the alignment
///
template <typename ScoreType>
struct GlobalAligner : public AlignerBase<ScoreType>
{
    GlobalAligner(
        const AlignmentScores<ScoreType>& scores) :
        AlignerBase<ScoreType>(scores)
    {}

    /// returns alignment path of query to reference
    template <typename SymIter>
    void
    align(
        const SymIter queryBegin, const SymIter queryEnd,
        const SymIter refBegin, const SymIter refEnd,
        AlignmentResult<ScoreType>& result) const;

private:

    // insert and delete are for query wrt reference
    struct ScoreVal
    {
        ScoreType match;
        ScoreType ins;
        ScoreType del;
    };

    struct PtrVal
    {
        uint8_t
        get(const AlignState::index_t i)
        {
            switch (i)
            {
            case AlignState::MATCH:
                return match;
            case AlignState::INSERT:
                return ins;
            case AlignState::DELETE:
                return del;
            default:
                assert(false && "Unexpected Index Value");
                return 0;
            }
        }

        uint8_t match : 2;
        uint8_t ins : 2;
        uint8_t del : 2;
    };

    // add the matrices here to reduce allocations over many alignment calls:
    typedef std::vector<ScoreVal> ScoreVec;
    mutable ScoreVec _score1;
    mutable ScoreVec _score2;
    mutable basic_matrix<PtrVal> _ptrMat;
};


#include "alignment/GlobalAlignerImpl.hh"
