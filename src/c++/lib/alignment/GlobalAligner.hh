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

#pragma once

#include "Alignment.hh"
#include "AlignmentScores.hh"

#include "blt_util/basic_matrix.hh"


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



struct AlignState
{
    enum index_t
    {
        MATCH,
        DELETE,
        INSERT,
        SIZE
    };

    static
    const char*
    label(const index_t i)
    {
        switch (i)
        {
        case MATCH:
            return "MATCH";
        case DELETE:
            return "DELETE";
        case INSERT:
            return "INSERT";
        default:
            return "UNKNOWN";
        }
    }
};



/// \brief Implementation of global alignment with affine gap costs
///
/// alignment outputs start positions and CIGAR-style alignment
/// expression of query to reference. Alignment of
/// seq1 is global -- seq1 can "fall-off" either end of the reference,
/// in this case, each unaligned position is scored as a mismatch and
/// the base is soft-clipped in the alignment
///
template <typename ScoreType>
struct GlobalAligner
{
    GlobalAligner(
        const AlignmentScores<ScoreType>& scores) :
        _scores(scores)
    {}

    /// returns alignment path of query to reference
    template <typename SymIter>
    void
    align(
        const SymIter queryBegin, const SymIter queryEnd,
        const SymIter refBegin, const SymIter refEnd,
        AlignmentResult<ScoreType>& result);

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


    const AlignmentScores<ScoreType> _scores;


    // insert and delete are for seq1 wrt seq2
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
                assert(! "Unexpected Index Value");
                return 0;
            }
        }

        uint8_t match : 2;
        uint8_t ins : 2;
        uint8_t del : 2;
    };

    // add the matrices here to reduce allocations over many alignment calls:
    typedef std::vector<ScoreVal> ScoreVec;
    ScoreVec _score1;
    ScoreVec _score2;
    basic_matrix<PtrVal> _ptrMat;
};


#include "alignment/GlobalAlignerImpl.hh"

