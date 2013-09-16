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
// <https://github.com/sequencing/licenses/>
//

///
/// \author Ole Schulz-Trieglaff
///

#pragma once

#include "alignment/Alignment.hh"
#include "alignment/GlobalJumpAligner.hh"

#include <cmath>

#include <boost/utility.hpp>

//#define DEBUG_RS


const double invlog10(1./(std::log(10.)));
const int phredScoreOffset = 33;

inline
pos_t
alignEnd(const Alignment& align)
{
    return (align.beginPos + ALIGNPATH::apath_ref_length(align.apath));
}


/// get begin position of alignment, accounting for a possibly reversed reference/alignment:
///
inline
pos_t
getAlignBeginOffset(
    const Alignment& align,
    const unsigned refSize,
    const bool isReversed)
{
    if (isReversed)
    {
        return (refSize - alignEnd(align));
    }
    else
    {
        return align.beginPos;
    }
}


/// get end position of alignment, accounting for a possibly reversed reference/alignment:
///
inline
pos_t
getAlignEndOffset(
    const Alignment& align,
    const unsigned refSize,
    const bool isReversed)
{
    if (isReversed)
    {
        return (refSize - align.beginPos);
    }
    else
    {
        return alignEnd(align);
    }
}


/// tests if prefix of aligned sequence matches target, returns length of alignment (zero if no match)
unsigned
hasAlignedPrefix(const Alignment& al, const unsigned minAlignContext = 0);


/// tests if suffix of aligned sequence matches target, returns length of alignment (zero if no match)
unsigned
hasAlignedSuffix(const Alignment& al, const unsigned minAlignContext = 0);


bool
bothEndsAligned(const Alignment& al, const unsigned minAlignContext = 0);


/// check a jump alignment for consistency (only one end aligning)
/// FIXME: not used, need to think what makes an alignment consistent
/// (how about : total number of matches shouldn't exceed sequence length?)
bool
isConsistentAlignment(const JumpAlignmentResult<int>& res, const unsigned minAlignContext = 0);


/// given a jump alignment and query sequence, return the bp1,insert and bp2 query sequences
///
/// the insert sequence is converted to fwd-strand by assuming it is "attached" to bp1
/// note this is targeted for debug-code only
///
void
getFwdStrandQuerySegments(
    const JumpAlignmentResult<int>& align,
    const std::string& querySeq,
    const bool isBp2AlignedFirst,
    const bool isBp1Reversed,
    const bool isBp2Reversed,
    std::string& bp1Seq,
    std::string& bp2Seq,
    std::string& insertSeq);


/// given a jump alignment and query sequence, return the insert sequence
///
/// the insert sequence is converted to fwd-strand by assuming it is "attached" to bp1
///
void
getFwdStrandInsertSegment(
    const JumpAlignmentResult<int>& align,
    const std::string& querySeq,
    const bool isBp1Reversed,
    std::string& insertSeq);


int
estimateBreakPointPos(const Alignment& al, const unsigned refOffset);



struct ReadScorer : private boost::noncopyable
{

    /** Instance getter
     *
    */
    static const ReadScorer& get()
    {
        static const ReadScorer rs;
        return rs;
    }

    double
    getSemiAlignedMetric(const ALIGNPATH::path_t& apath, const uint8_t* qual) const;

private:
    explicit
    ReadScorer();
    ~ReadScorer() {}

    enum { MAX_Q = 128 };
    const int _qmin;
    double _logpcorrectratio[MAX_Q];
};
