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

#include "alignment/AlignmentUtil.hh"
#include "blt_util/seq_util.hh"

#include <cassert>

//#define DEBUG_RS



#ifdef DEBUG_RS
#include "blt_util/log.hh"
#include <iostream>
#endif



/// tests if prefix of aligned sequence matches target, returns length of alignment (zero if no match)
static
unsigned
hasAlignedPrefix(const Alignment& al, const unsigned minAlignContext)
{
    if (al.apath.empty()) return false;
    unsigned alignLen(0);
    if (al.apath[0].type == ALIGNPATH::MATCH && al.apath[0].length >= minAlignContext)
    {
        alignLen = al.apath[0].length;
    }
    return alignLen;
}



// tests if suffix of aligned sequence matches target, returns length of alignment (zero if no match)
static
unsigned
hasAlignedSuffix(const Alignment& al, const unsigned minAlignContext)
{
    if (al.apath.empty()) return false;
    size_t apLen = al.apath.size();
    unsigned alignLen(0);
    if (al.apath[apLen-1].type == ALIGNPATH::MATCH && al.apath[apLen-1].length >= minAlignContext)
    {
        alignLen = al.apath[apLen-1].length;
    }
    return alignLen;
}


#if 0
bool
bothEndsAligned(const Alignment& al, const unsigned minAlignContext)
{
    return (hasAlignedPrefix(al,minAlignContext) && hasAlignedSuffix(al,minAlignContext));
}
#endif



// check a jump alignment for consistency (only one end aligning)
// FIXME: not used, need to think what makes an alignment consistent
// (how about : total number of matches shouldn't exceed sequence length?)
//bool
//isConsistentAlignment(const JumpAlignmentResult<int>& res, const unsigned /*minAlignContext = 0*/)
//{
//    // not consistent if both unaligned
//    if (! (res.align1.isAligned() && res.align2.isAligned()) ) return false;
//
//    return true;
//}



void
getFwdStrandQuerySegments(
    const JumpAlignmentResult<int>& align,
    const std::string& querySeq,
    const bool isBp2AlignedFirst,
    const bool isBp1Reversed,
    const bool isBp2Reversed,
    std::string& bp1Seq,
    std::string& bp2Seq,
    std::string& insertSeq)
{
    const unsigned align1Size(apath_read_length(align.align1.apath));
    const unsigned insertSize(align.jumpInsertSize);
    const unsigned insertOffset(align1Size + insertSize);
    const unsigned align2Size(apath_read_length(align.align2.apath));
    const unsigned align2Offset(insertOffset + align2Size);

    assert(querySeq.size() == align2Offset);

    bp1Seq = querySeq.substr(0,align1Size);
    insertSeq = querySeq.substr(align1Size,insertSize);
    bp2Seq = querySeq.substr(insertOffset,align2Size);

    if (isBp2AlignedFirst)
    {
        std::swap(bp1Seq,bp2Seq);
    }

    if (isBp1Reversed)
    {
        reverseCompStr(bp1Seq);
        reverseCompStr(insertSeq);
    }

    if (isBp2Reversed)
    {
        reverseCompStr(bp2Seq);
    }
}



void
getExtendedContig(
    const JumpAlignmentResult<int>& align,
    const std::string& querySeq,
    const std::string* ref1Seq,
    const std::string* ref2Seq,
    std::string& extendedContig)
{
    // extend aligned ref1Seq to the beginning
    const std::string align1RefSeqExt = ref1Seq->substr(0, align.align1.beginPos);
    // extend aligned ref2Seq to the end
    const unsigned align2RefEnd = align.align2.beginPos + apath_ref_length(align.align2.apath);
    const std::string align2RefSeqExt = ref2Seq->substr(align2RefEnd, (ref2Seq->size()-align2RefEnd));

    extendedContig = align1RefSeqExt + querySeq + align2RefSeqExt;
}



void
getExtendedContig(
    const AlignmentResult<int>& alignment,
    const std::string& querySeq,
    const std::string* refSeq,
    std::string& extendedContig)
{
    // extend aligned refSeq to the beginning
    const std::string alignRefSeqExt1 = refSeq->substr(0, alignment.align.beginPos);
    // extend aligned refSeq to the end
    const unsigned alignRefEnd = alignment.align.beginPos + apath_ref_length(alignment.align.apath);
    const std::string alignRefSeqExt2 = refSeq->substr(alignRefEnd, (refSeq->size()-alignRefEnd));

    extendedContig = alignRefSeqExt1 + querySeq + alignRefSeqExt2;
}



void
getFwdStrandInsertSegment(
    const JumpAlignmentResult<int>& align,
    const std::string& querySeq,
    const bool isBp1Reversed,
    std::string& insertSeq)
{
    const unsigned align1Size(apath_read_length(align.align1.apath));
    const unsigned insertSize(align.jumpInsertSize);

    insertSeq = querySeq.substr(align1Size,insertSize);

    if (isBp1Reversed) reverseCompStr(insertSeq);
}



int
estimateBreakPointPos(const Alignment& al, const unsigned refOffset)
{
    // -1 means no breakpoint found
    int breakPointPosEstimate(-1);

    const unsigned prefAlLen = hasAlignedPrefix(al, 0);
    const unsigned suffAlLen = hasAlignedSuffix(al, 0);

    if (! (prefAlLen || suffAlLen) )
    {
        return breakPointPosEstimate;
    }

    if (prefAlLen)
    {
        breakPointPosEstimate = refOffset + al.beginPos /*+ prefAlLen*/;
    }
    else if (suffAlLen)
    {
        breakPointPosEstimate = refOffset + alignEnd(al) /*- suffAlLen*/;
    }

    assert(breakPointPosEstimate>0);

    return breakPointPosEstimate;
}
