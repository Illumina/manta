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

///
/// \author Chris Saunders
///

#pragma once

#include "JumpAlignerBase.hh"



/// \brief a method to align a contig to two references
///
/// the alignment can make a single jump from reference1 to reference2
///
/// transition from insert to delete is free and allowed, but not reverse
/// transition from/to jump to insert is free and allowed TODO: more restrictive
///
template <typename ScoreType>
struct GlobalJumpAligner : public JumpAlignerBase<ScoreType>
{
    GlobalJumpAligner(
        const AlignmentScores<ScoreType>& scores,
        const ScoreType jumpScore) :
        JumpAlignerBase<ScoreType>(scores,jumpScore)
    {}

    /// returns alignment path of query to reference
    template <typename SymIter>
    void
    align(
        const SymIter queryBegin, const SymIter queryEnd,
        const SymIter ref1Begin, const SymIter ref1End,
        const SymIter ref2Begin, const SymIter ref2End,
        JumpAlignmentResult<ScoreType>& result) const;

private:

    // insert and delete are for seq1 wrt seq2
    struct ScoreVal
    {
        ScoreType match;
        ScoreType ins;
        ScoreType del;
        ScoreType jump;
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
            case AlignState::JUMP:
                return jump;
            default:
                assert(false && "Unexpected Index Value");
                return 0;
            }
        }

        /// pack 2x4 bits into 1 byte:
        uint8_t match : 2;
        uint8_t ins : 2;
        uint8_t del : 2;
        uint8_t jump : 2;
    };

    // add the matrices here to reduce allocations over many alignment calls:
    typedef std::vector<ScoreVal> ScoreVec;
    mutable ScoreVec _score1;
    mutable ScoreVec _score2;

    typedef basic_matrix<PtrVal> PtrMat;
    mutable PtrMat _ptrMat1;
    mutable PtrMat _ptrMat2;
};


#include "alignment/GlobalJumpAlignerImpl.hh"
