
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
      _data(rowCount*colCount)
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
    empty() { return _data.empty(); }

    size_t
    size() { return _data.size(); }

    iterator
    begin() { return _data.begin(); }

    const_iterator
    begin() const { return _data.begin(); }

    iterator
    end() { return _data.end(); }

    const_iterator
    end() const { return _data.end(); }

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


struct AlignmentResult {

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
        SymIter begin1, SymIter end1,
        SymIter begin2, SymIter end2,
        AlignmentResult& result);

private:

    void
    max3(
        ScoreType& max,
        uint8_t& which,
        const ScoreType v0,
        const ScoreType v1,
        const ScoreType v2)
    {
        max=v0;
        which=0;
        if (v1>v0)
        {
            max=v1;
            which=1;
        }
        if (v2>max)
        {
            max=v2;
            which=2;
        }
    }


    const AlignmentScores<ScoreType> _scores;

    // insert and delete are for seq1 wrt seq2
    struct ScoreVal
    {
        ScoreType match;
        ScoreType ins;
        ScoreType del;
    };

    struct WhichVal
    {
        uint8_t match : 2;
        uint8_t ins : 2;
        uint8_t del : 2;
    };

    BasicMatrix<ScoreVal> _scoreMat;
    BasicMatrix<WhichVal> _whichMat;
};


#include "alignment/GlobalAlignerImpl.hh"

