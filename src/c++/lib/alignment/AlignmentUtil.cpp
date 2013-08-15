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
/// \author Ole Schulz-Trieglaff
///

#include "alignment/AlignmentUtil.hh"

#include <cassert>


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
bool
isConsistentAlignment(const JumpAlignmentResult<int>& res, const unsigned /*minAlignContext = 0*/)
{
    // not consistent if both unaligned
    if (! (res.align1.isAligned() && res.align2.isAligned()) ) return false;

    return true;
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

    if (prefAlLen) {
        breakPointPosEstimate = refOffset + al.alignStart + prefAlLen;
    } else if (suffAlLen) {
        breakPointPosEstimate = refOffset + alignEnd(al) - suffAlLen;
    }

    assert(breakPointPosEstimate>0);

    return breakPointPosEstimate;
}
