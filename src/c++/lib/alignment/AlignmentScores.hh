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

#pragma once


template <typename ScoreType>
struct AlignmentScores
{
    static const ScoreType badVal = -10000;

    AlignmentScores(
        ScoreType initMatch,
        ScoreType initMismatch,
        ScoreType initOpen,
        ScoreType initExtend,
        ScoreType initIntronOpen,
        ScoreType initOffEdge,
        ScoreType initIntronOffEdge) :
        match(initMatch),
        mismatch(initMismatch),
        open(initOpen),
        extend(initExtend),
        intronOpen(initIntronOpen),
        offEdge(initOffEdge),
        intronOffEdge(initIntronOffEdge)
    {}

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
        intronOpen(badVal),
        offEdge(initOffEdge),
        intronOffEdge(badVal)
    {}

    const ScoreType match; ///< match score
    const ScoreType mismatch; ///< mismatch score (should be negative)
    const ScoreType open; ///< gap open, gap of length N is scored (open + N * extend) (should be negative)
    const ScoreType extend; ///< gap extend, gap of length N is scored (open + N * extend) (should be negative or zero)
    const ScoreType intronOpen; ///< gap open for introns (i.e. deletions starting with splice motif) (should be negative)
    const ScoreType offEdge; ///< score applied when query goes off the end of an edge (should be negative)
    const ScoreType intronOffEdge; ///As offEdge but only of the last aligned bases match a splice motif (should be negative) todo: not implemented
};
