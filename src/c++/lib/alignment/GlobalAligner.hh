// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Manta
// Copyright (c) 2013-2015 Illumina, Inc.
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

#include "SingleRefAlignerShared.hh"


/// \brief Implementation of global alignment with affine gap costs
///
/// alignment outputs start positions and CIGAR-style alignment
/// expression of query to reference. Alignment of
/// query is global -- query can "fall-off" either end of the reference,
/// in this case, each unaligned position is given an "off-edge" score and
/// the base is soft-clipped in the alignment
///
/// transition from insert to delete is free and allowed, but not reverse
///
template <typename ScoreType>
struct GlobalAligner : public SingleRefAlignerBase<ScoreType>
{
    GlobalAligner(
        const AlignmentScores<ScoreType>& scores) :
        SingleRefAlignerBase<ScoreType>(scores)
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
        get(const AlignState::index_t i) const
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
