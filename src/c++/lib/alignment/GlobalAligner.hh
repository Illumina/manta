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

/// based off of ELAND implementation by Tony Cox

#pragma once

#include "blt_util/align_path.hh"

#include <vector>



/// very simple matrix implementation, row major
template <typename T>
struct BasicMatrix
{
    typedef typename std::vector<T> data_t;
    typedef typename data_t::iterator iterator;
    typedef typename data_t::const_iterator const_iterator;

    BasicMatrix(
        const unsigned rowCount = 0,
        const unsigned colCount = 0) :
        _colCount(colCount),
        _data(rowCount* colCount)
    {}

    void
    resize(
        const unsigned rowCount,
        const unsigned colCount)
    {
        _colCount=colCount;
        _data.resize(rowCount*colCount);
    }

    T&
    val(const unsigned row,
        const unsigned col)
    {
        return _data[(row*_colCount+col)];
    }

    const T&
    val(const unsigned row,
        const unsigned col) const
    {
        return _data[(row*_colCount+col)];
    }

    bool
    empty()
    {
        return _data.empty();
    }

    size_t
    size()
    {
        return _data.size();
    }

    iterator
    begin()
    {
        return _data.begin();
    }

    const_iterator
    begin() const
    {
        return _data.begin();
    }

    iterator
    end()
    {
        return _data.end();
    }

    const_iterator
    end() const
    {
        return _data.end();
    }

private:
    unsigned _colCount;
    std::vector<T> _data;
};



template <typename ScoreType>
struct AlignmentScores
{
    AlignmentScores(
        ScoreType initMatch,
        ScoreType initMismatch,
        ScoreType initOpen,
        ScoreType initExtend) :
        match(initMatch),
        mismatch(initMismatch),
        open(initOpen),
        extend(initExtend)
    {}

    const ScoreType match;
    const ScoreType mismatch;
    const ScoreType open;
    const ScoreType extend;
};


template <typename ScoreType>
struct AlignmentResult
{

    void
    clear()
    {
        score=0;
        alignStart=0;
        apath.clear();
    }

    ScoreType score;
    unsigned alignStart;
    ALIGNPATH::path_t apath;
};



/**
 * @class Aligner
 *
 * @brief Implementation of global alignment with affine gap costs.
 *
 */
template <typename ScoreType>
struct GlobalAligner
{
    GlobalAligner(AlignmentScores<ScoreType>& scores) :
        _scores(scores)
    {}

    /// returns alignment path of 1 to 2
    template <typename SymIter>
    void
    align(
        const SymIter begin1, const SymIter end1,
        const SymIter begin2, const SymIter end2,
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

    void
    updatePath(ALIGNPATH::path_t& path,
               ALIGNPATH::path_segment& ps,
               ALIGNPATH::align_t atype)
    {
        if (ps.type == atype) return;
        if (ps.type != ALIGNPATH::NONE) path.push_back(ps);
        ps.type = atype;
        ps.length = 0;
    }


    const AlignmentScores<ScoreType> _scores;

    enum matrix_t
    {
        MATCH,
        INSERT,
        DELETE,
        SIZE
    };

    // insert and delete are for seq1 wrt seq2
    struct ScoreVal
    {
        ScoreType
        get(const matrix_t i)
        {
            switch (i)
            {
            case MATCH:
                return match;
            case INSERT:
                return ins;
            case DELETE:
                return del;
            }
        }

        ScoreType match;
        ScoreType ins;
        ScoreType del;
    };

    struct PtrVal
    {
        uint8_t
        get(const matrix_t i)
        {
            switch (i)
            {
            case MATCH:
                return match;
            case INSERT:
                return ins;
            case DELETE:
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

    BasicMatrix<ScoreVal> _scoreMat;
    BasicMatrix<PtrVal> _ptrMat;
};


#include "alignment/GlobalAlignerImpl.hh"

