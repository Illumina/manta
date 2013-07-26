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

#pragma once


template <typename ScoreType>
struct AlignmentScores
{
    AlignmentScores(
        ScoreType initMatch,
        ScoreType initMismatch,
        ScoreType initOpen,
        ScoreType initExtend,
        ScoreType initOffEdge) :
        match(initMatch),
        mismatch(initMismatch),
        open(initOpen),
        extend(initExtend),
        offEdge(initOffEdge)
    {}

    const ScoreType match;
    const ScoreType mismatch;
    const ScoreType open; /// gap open, gap of length N is scored (open + N * extend)
    const ScoreType extend; /// gap extend, gap of length N is scored (open + N * extend)
    const ScoreType offEdge; /// score applied when query goes off the end of an edge
};
