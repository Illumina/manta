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
#include <math.h>


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
        breakPointPosEstimate = refOffset + al.beginPos /*+ prefAlLen*/;
    } else if (suffAlLen) {
        breakPointPosEstimate = refOffset + alignEnd(al) /*- suffAlLen*/;
    }

    assert(breakPointPosEstimate>0);

    return breakPointPosEstimate;
}




/// returns log(1+x), switches to special libc function when abs(x) is small
///
static
double
log1p_switch(const double x) {
    // better number??
    static const double smallx_thresh(0.01);
    if (std::abs(x)<smallx_thresh) {
        return log1p(x);
    } else {
        return std::log(1+x);
    }
}

static
double convertPhredToProbError(int qv) {
    return std::min(1.,std::pow(10.,(-static_cast<double>(qv)/10.)));
}



ReadScorer::
ReadScorer()
    : _qmin(phredScoreOffset) {

#ifdef DEBUG_SU
    std::cout << "Filling logpcorrectratio table" << std::endl;
#endif
    _logpcorrectratio[_qmin] = 0;
    for (int i(_qmin+1); i<MAX_Q; ++i) {
        const double eprob(convertPhredToProbError(i-phredScoreOffset));
#ifdef DEBUG_SU
        std::cout << "i=" << i << " " << log1p_switch(-eprob) << " " << std::log(eprob/3.) << std::endl;
#endif
        _logpcorrectratio[i] = log1p_switch(-eprob) - std::log(eprob/3.);
    }

#ifdef DEBUG_SU
    std::cout << "Readscorer dumping _logpcorrectratio : " << std::endl;
    for (int i(_qmin); i<MAX_Q; ++i) {
        std::cout << "_logpcorrectratio[" <<  i << "] = " << _logpcorrectratio[i] << std::endl;
    }
#endif
}



double
ReadScorer::
getSemiAlignedMetric(const std::string& matchDescriptor, const std::string& qualities) {

    int posInRead = 0;
    int tmp = 0;
    double alignScore(0.);
    bool inGap(false);

/*    if (!isValidQualString(qualityString)) {
#ifdef DEBUG_SU
        std::cout << "!isValidQualString(qualityString)" << std::endl;
#endif
        return alignScore;
    }
*/

#ifdef DEBUG_SU
    std::cout << "getAlignmentScore type=" << matchDescriptor << std::endl;
    std::cout << "getAlignmentScore qual=" << qualities << std::endl;
#endif
    const unsigned nt(matchDescriptor.size());
    for (unsigned i=0; i<nt; ++i)
    {
#ifdef DEBUG_SU
        std::cout << "getAlignmentScore : i = " << i << " alignScore = " << alignScore << std::endl;
#endif
        // ASCII: 48 = 0, 57=9
        if ((int) matchDescriptor[i] <= 57 && (int) matchDescriptor[i] >= 48)
        {
#ifdef DEBUG_SU
            std::cout << "getAlignmentScore : i = " << i << " alignScore = " << alignScore << std::endl;
#endif
            tmp *= 10;
            tmp += (matchDescriptor[i] - '0');
        }
        else if (matchDescriptor[i] == 'A' || matchDescriptor[i] == 'C' || matchDescriptor[i] == 'G' || matchDescriptor[i] == 'T')
        {
            posInRead += tmp;
            tmp=0;
            // Bases in gap -> dels : ignoring for scoring purposes
            //                        & do not advance the read position.
            if (!inGap) {
                // mismatch
#ifdef DEBUG_SU
                std::cerr << "getAlignmentScore: " << posInRead << " " << _logpcorrectratio[static_cast<int>(qualities[posInRead])]
                          << " " << qualities[posInRead] << " " << static_cast<int>(qualities[posInRead]) << std::endl;
#endif
                alignScore +=  _logpcorrectratio[static_cast<int>(qualities[posInRead])];
                ++posInRead;
            }
        } else if (matchDescriptor[i] == 'N') {
            posInRead += tmp;
            tmp=0;
            // what should N in the reference be?? for now it counts as a match
            ++posInRead;
        } else if (matchDescriptor[i]=='^') {
            posInRead += tmp;
            tmp=0;
            inGap = true;
        } else if (matchDescriptor[i] == '$') {
            // If tmp is non-zero here, it indicates an insertion.
            // FIXME : probably should update score using open & extend penalties.
            // Currently just advancing the read position.
            posInRead += tmp;
            tmp=0;
            inGap = false;
        } else {
        	// this
            assert(!"ERROR: Unexpected match descriptor!");
        }
    }
    return alignScore;


}
