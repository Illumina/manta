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

///
/// \author Chris Saunders and Xiaoyu Chen
///

#include "SVScorePairRefProcessor.hh"
#include "SVScorerShared.hh"
#include "blt_util/bam_record_util.hh"
#include "common/Exceptions.hh"
#include "manta/SVCandidateUtil.hh"

#include <cassert>

#include <sstream>


/// standard debug output for this file:
//#define DEBUG_PAIR

/// ridiculous debug output for this file:
//#define DEBUG_MEGAPAIR

#ifdef DEBUG_PAIR
#include "blt_util/log.hh"
#endif



void
SVScorePairRefProcessor::
processClearedRecord(
    const bam_record& bamRead)
{
    using namespace illumina::common;

    assert(bamParams.isSet);

    const pos_t refPos(bamRead.pos()-1);
    if (! bamParams.interval.range.is_pos_intersect(refPos)) return;

#ifdef DEBUG_MEGAPAIR
    log_os << __FUNCTION__ << ": read: " << bamRead << "\n";
#endif

    /// check if fragment is too big or too small:
    const int templateSize(std::abs(bamRead.template_size()));
    if (templateSize < bamParams.minFrag) return;
    if (templateSize > bamParams.maxFrag) return;

    // count only from the down stream reads
    const bool isFirstBamRead(isFirstRead(bamRead));

    // get fragment range:
    pos_t fragBeginRefPos(refPos);
    if (! isFirstBamRead)
    {
        fragBeginRefPos=bamRead.mate_pos()-1;
    }

    const pos_t fragEndRefPos(fragBeginRefPos+templateSize);

    if (fragBeginRefPos > fragEndRefPos)
    {
        std::ostringstream oss;
        oss << "ERROR: Failed to parse fragment range from bam record. Frag begin,end: " << fragBeginRefPos << " " << fragEndRefPos << " bamRecord: " << bamRead << "\n";
        BOOST_THROW_EXCEPTION(LogicException(oss.str()));
    }

    {
        const pos_t fragOverlap(std::min((1+svParams.centerPos-fragBeginRefPos), (fragEndRefPos-svParams.centerPos)));
#ifdef DEBUG_MEGAPAIR
        log_os << __FUNCTION__ << ": frag begin/end/overlap: " << fragBeginRefPos << " " << fragEndRefPos << " " << fragOverlap << "\n";
#endif
        if (fragOverlap < pairOpt.minFragSupport) return;
    }

    SVFragmentEvidence& fragment(evidence.getSample(bamParams.isTumor)[bamRead.qname()]);

    static const bool isShadow(false);

    SVFragmentEvidenceRead& evRead(fragment.getRead(bamRead.is_first()));
    setReadEvidence(svParams.minMapQ, svParams.minTier2MapQ, bamRead, isShadow, evRead);

    setAlleleFrag(*bamParams.fragDistroPtr, templateSize, fragment.ref.getBp(isBp1));
}
