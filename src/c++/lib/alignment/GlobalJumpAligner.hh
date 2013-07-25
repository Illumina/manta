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

/// \author Chris Saunders

#pragma once

#include "Alignment.hh"
#include "AlignmentScores.hh"

#include "blt_util/basic_matrix.hh"

#include <cassert>


/// represents alignment of a query sequence which can switch over from reference1 to refernece2
///
/// an empty alignment to one reference indicates that the entire alignment is the other reference
///
template <typename ScoreType>
struct JumpAlignmentResult
{
    JumpAlignmentResult()
    {
        clear();
    }

    void
    clear()
    {
        score=0;
        jumpInsertSize=0;
        align1.clear();
        align2.clear();
    }


    ScoreType score;
    unsigned jumpInsertSize;
    Alignment align1;
    Alignment align2;
};




struct AlignState
{
    enum index_t
    {
        MATCH,
        DELETE,
        INSERT,
        JUMP,
        SIZE
    };

    static
    const char*
    label(const index_t i)
    {
        switch(i)
        {
        case MATCH:  return "MATCH";
        case DELETE: return "DELETE";
        case INSERT: return "INSERT";
        case JUMP: return "JUMP";
        default:     return "UNKNOWN";
        }
    }
};



/// \brief a method to align a contig to two references, where the alignment can make
/// a single jump from reference one to reference 2
///
template <typename ScoreType>
struct GlobalJumpAligner
{
    GlobalJumpAligner(
        const AlignmentScores<ScoreType>& scores,
        const ScoreType jumpScore) :
        _scores(scores),
        _jumpScore(jumpScore)
    {}

    /// returns alignment path of query to reference
    template <typename SymIter>
    void
    align(
        const SymIter queryBegin, const SymIter queryEnd,
        const SymIter ref1Begin, const SymIter ref1End,
        const SymIter ref2Begin, const SymIter ref2End,
        JumpAlignmentResult<ScoreType>& result);

private:

    uint8_t
    max3(
        ScoreType& max,
        const ScoreType v0,
        const ScoreType v1,
        const ScoreType v2)
    {
        max=v0;
        uint8_t ptr=0;
        if (v1>v0)
        {
            max=v1;
            ptr=1;
        }
        if (v2>max)
        {
            max=v2;
            ptr=2;
        }
        return ptr;
    }

    uint8_t
    max4(
        ScoreType& max,
        const ScoreType v0,
        const ScoreType v1,
        const ScoreType v2,
        const ScoreType v3)
    {
        max=v0;
        uint8_t ptr=0;
        if (v1>v0)
        {
            max=v1;
            ptr=1;
        }
        if (v2>max)
        {
            max=v2;
            ptr=2;
        }
        if (v3>max)
        {
            max=v3;
            ptr=3;
        }
        return ptr;
    }

    const AlignmentScores<ScoreType> _scores;
    const ScoreType _jumpScore;

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
            case AlignState::JUMP:
                return jump;
            default:
                assert(! "Unexpected Index Value");
                return 0;
            }
        }

        uint8_t match : 2;
        uint8_t ins : 2;
        uint8_t del : 2;
        uint8_t jump : 2;
    };

    // add the matrices here to reduce allocations over many alignment calls:
    typedef std::vector<ScoreVal> ScoreVec;
    ScoreVec _score1;
    ScoreVec _score2;

    typedef basic_matrix<PtrVal> PtrMat;
    PtrMat _ptrMat1;
    PtrMat _ptrMat2;
};


#include "alignment/GlobalJumpAlignerImpl.hh"

