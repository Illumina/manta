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

#define DEBUG_RS

#ifdef DEBUG_RS
#include "blt_util/log.hh"
//#include <iostream>
#endif

// tests if prefix of aligned sequence matches target, returns length of alignment (zero if no match)
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


bool
bothEndsAligned(const Alignment& al, const unsigned minAlignContext)
{
    return (hasAlignedPrefix(al,minAlignContext) && hasAlignedSuffix(al,minAlignContext));
}



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

    unsigned prefAlLen = hasAlignedPrefix(al);
    unsigned suffAlLen = hasAlignedSuffix(al);

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




/// returns log(1+x), switches to special libc function when abs(x) is small
///
static
double
log1p_switch(const double x)
{
    // better number??
    static const double smallx_thresh(0.01);
    if (std::abs(x)<smallx_thresh)
    {
        return log1p(x);
    }
    else
    {
        return std::log(1+x);
    }
}

static
double convertPhredToProbError(int qv)
{
    return std::min(1.,std::pow(10.,(-static_cast<double>(qv)/10.)));
}



ReadScorer::
ReadScorer()
    : _qmin(phredScoreOffset)
{

#ifdef DEBUG_RS
    std::cout << "Filling logpcorrectratio table" << std::endl;
#endif
    _logpcorrectratio[_qmin] = 0;
    for (int i(_qmin+1); i<MAX_Q; ++i)
    {
        const double eprob(convertPhredToProbError(i-phredScoreOffset));
#ifdef DEBUG_RS
        std::cout << "i=" << i << " " << log1p_switch(-eprob) << " " << std::log(eprob/3.) << std::endl;
#endif
        _logpcorrectratio[i] = log1p_switch(-eprob) - std::log(eprob/3.);
    }

#ifdef DEBUG_RS
    std::cout << "Readscorer dumping _logpcorrectratio : " << std::endl;
    for (int i(_qmin); i<MAX_Q; ++i)
    {
        std::cout << "_logpcorrectratio[" <<  i << "] = " << _logpcorrectratio[i] << std::endl;
    }
#endif
}


double
ReadScorer::
getSemiAlignedMetric(
    const ALIGNPATH::path_t& apath,
    const uint8_t* qual) const
{
    using namespace ALIGNPATH;

    unsigned posInRead = 0;
    double score(0.);

#ifdef DEBUG_RS
    log_os << "getAlignmentScore apath=" << apath << "\n";
#endif

    const unsigned nt(apath.size());
    for (unsigned i=0; i<nt; ++i)
    {
        const path_segment& ps(apath[i]);
#ifdef DEBUG_RS
        log_os << "getSemiAlignedMetric : i=" << i << " posInRead=" << posInRead << " score=" << score << "\n";
#endif
        if((ps.type==SOFT_CLIP) || (ps.type==SEQ_MISMATCH))
        {
            for(unsigned j(0);j<ps.length;++j)
            {
#ifdef DEBUG_RS
                log_os << "getAlignmentScore: " << (posInRead+j) << " " << _logpcorrectratio[qual[posInRead+j]]
                       << " " << qual[posInRead] << "\n";
#endif
                score +=  _logpcorrectratio[qual[posInRead+j]];
            }
         }
        if(is_segment_type_read_length(ps.type)) posInRead += ps.length;
   } // for
   return score;
}


