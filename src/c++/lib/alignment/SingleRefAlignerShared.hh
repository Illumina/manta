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

/// common boilerplate for single-reference sequence aligners

#pragma once

#include "AlignerBase.hh"
#include "AlignerUtil.hh"
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



template <typename ScoreType>
struct SingleRefAlignerBase : public AlignerBase<ScoreType>
{
    SingleRefAlignerBase(
        const AlignmentScores<ScoreType>& scores) :
        AlignerBase<ScoreType>(scores)
    {}

protected:

    /// returns alignment path of query to reference
    template <typename SymIter, typename MatrixType>
    void
    backTraceAlignment(
        const SymIter queryBegin, const SymIter queryEnd,
        const SymIter refBegin, const SymIter refEnd,
        const size_t querySize, const size_t refSize,
        const MatrixType& ptrMatrix,
        const BackTrace<ScoreType>& btraceInput,
        AlignmentResult<ScoreType>& result) const;
};


#include "SingleRefAlignerSharedImpl.hh"
