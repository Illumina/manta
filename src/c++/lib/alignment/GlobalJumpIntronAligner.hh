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

///
/// \author Chris Saunders and Felix Schlesinger
///

#pragma once

#include "JumpAlignerBase.hh"


/// \brief a method to align a contig to two references
///
/// the alignment can make a single jump from reference1 to reference2
///
template <typename ScoreType>
struct GlobalJumpIntronAligner : public JumpAlignerBase<ScoreType>
{
    GlobalJumpIntronAligner(
        const AlignmentScores<ScoreType>& scores,
        const ScoreType jumpScore,
        const ScoreType intronOpenScore,
        const ScoreType intronOffEdgeScore) :
        JumpAlignerBase<ScoreType>(scores,jumpScore),
        _intronOpenScore(intronOpenScore),
        _intronOffEdgeScore(intronOffEdgeScore)
    {}

    /// returns alignment path of query to reference
    template <typename SymIter>
    void
    align(
        const SymIter queryBegin, const SymIter queryEnd,
        const SymIter ref1Begin, const SymIter ref1End,
        const SymIter ref2Begin, const SymIter ref2End,
        bool ref1Fw, bool ref2Fw,
        JumpAlignmentResult<ScoreType>& result) const;

private:

    // insert and delete are for seq1 wrt seq2
    struct ScoreVal
    {
        ScoreType match;
        ScoreType ins;
        ScoreType del;
        ScoreType jump;
        ScoreType intron;
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
            case AlignState::SPLICE:
                return intron;
            default:
                assert(false && "Unexpected Index Value");
                return 0;
            }
        }

        /// pack 3x5 bits into a single uint16_t:
        uint16_t match : 3;
        uint16_t ins : 3;
        uint16_t del : 3;
        uint16_t jump : 3;
        uint16_t intron : 3;
    };

    const ScoreType _intronOpenScore; ///< gap open for introns (i.e. deletions starting with splice motif) (should be negative)
    const ScoreType _intronOffEdgeScore; ///< As offEdge but only of the last aligned bases match a splice motif (should be negative) todo: not implemented

    // add the matrices here to reduce allocations over many alignment calls:
    typedef std::vector<ScoreVal> ScoreVec;
    mutable ScoreVec _score1;
    mutable ScoreVec _score2;

    typedef basic_matrix<PtrVal> PtrMat;
    mutable PtrMat _ptrMat1;
    mutable PtrMat _ptrMat2;
};


#include "alignment/GlobalJumpIntronAlignerImpl.hh"
