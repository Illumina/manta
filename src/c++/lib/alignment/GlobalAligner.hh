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
        switch(i)
        {
        case MATCH:  return "MATCH";
        case DELETE: return "DELETE";
        case INSERT: return "INSERT";
        default:     return "UNKNOWN";
        }
    }
};



/// \brief Implementation of global alignment with affine gap costs
///
/// alignment outputs start positions and CIGAR-style alignment
/// expression of seq1 (query) to seq2 (referenece). Alignment of
/// seq1 is global -- seq1 can "fall-off" either end of the reference,
/// in this case, each unaligned postiion is scored as a mismatch and
/// the base is soft-clipped in the alignment
///
template <typename ScoreType>
struct GlobalAligner
{
    GlobalAligner(AlignmentScores<ScoreType>& scores) :
        _scores(scores)
    {}

    /// returns alignment path of 1 (query) to 2 (reference)
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
    BasicMatrix<ScoreVal> _scoreMat;
    BasicMatrix<PtrVal> _ptrMat;
};


#include "alignment/GlobalAlignerImpl.hh"

