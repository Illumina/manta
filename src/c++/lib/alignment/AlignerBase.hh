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

///
/// \author Chris Saunders
///

#pragma once

#include "alignment/AlignmentScores.hh"
#include "blt_util/align_path.hh"

#include "boost/foreach.hpp"

#include <cassert>


/// base class containing a few shared aligner functions
///
template <typename ScoreType>
struct AlignerBase
{
    AlignerBase(
        const AlignmentScores<ScoreType>& scores) :
        _scores(scores)
    {}

    /// read-only access to the aligner's scores:
    const AlignmentScores<ScoreType>&
    getScores() const
    {
        return _scores;
    }

    /// recover a path alignment score without aligning, requires SEQ_MATCH style CIGAR
    ///
    ScoreType
    getPathScore(
        const ALIGNPATH::path_t& apath,
        const bool isScoreOffEdge = true) const
    {
        using namespace ALIGNPATH;

        ScoreType val(0);

        BOOST_FOREACH(const path_segment& ps, apath)
        {
            bool isIndel(false);
            switch (ps.type)
            {
            case MATCH:
                assert(! "Unexpected MATCH segment"); // if MATCH segments exist, then you're using the wrong type of CIGAR for this function
                break;
            case SEQ_MATCH:
                val += (_scores.match * ps.length);
                isIndel = false;
                break;
            case SEQ_MISMATCH:
                val += (_scores.mismatch * ps.length);
                isIndel = false;
                break;
            case INSERT:
            case DELETE:
                if (! isIndel) val += _scores.open;
                val += (_scores.extend * ps.length);
                isIndel = true;
                break;
            case SOFT_CLIP:
                if (isScoreOffEdge) val += (_scores.offEdge * ps.length);
                isIndel = false;
                break;
            default:
                break;
            }
        }
        return val;
    }

protected:

    static
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
};
