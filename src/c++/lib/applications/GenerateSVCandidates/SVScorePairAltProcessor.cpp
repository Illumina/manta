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

#include "SVScorePairAltProcessor.hh"
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
SVScorePairAltProcessor::
checkInput(
    const SVCandidate& sv)
{
    using namespace illumina::common;

    // this class is designed for simple alts only:
    assert(sv.bp1.interval.tid == sv.bp2.interval.tid);
    assert(getSVType(sv) == SV_TYPE::INDEL);

    /// In case of breakend microhomology approximate the breakend as a point event at the center of the possible range:
    const pos_t centerPos1 = (sv.bp1.interval.range.center_pos());
    const pos_t centerPos2 = (sv.bp2.interval.range.center_pos());
    if (centerPos2 <= centerPos1)
    {
        std::ostringstream oss;
        oss << "ERROR: Unexpected breakend orientation in pair support routine for sv: " << sv << "\n";
        BOOST_THROW_EXCEPTION(LogicException(oss.str()));
    }
}



void
SVScorePairAltProcessor::
processClearedRecord(
    const bam_record& bamRead)
{
    using namespace illumina::common;

    assert(bparams.isSet);

    const pos_t refPos(bamRead.pos()-1);
    if (! bparams.interval.range.is_pos_intersect(refPos)) return;

#ifdef DEBUG_MEGAPAIR
    log_os << __FUNCTION__ << ": read: " << bamRead << "\n";
#endif

    /// check if fragment is too big or too small:
    const int templateSize(std::abs(bamRead.template_size()));
    const int altTemplateSize(templateSize-iparams.altShift);
    if (altTemplateSize < bparams.minFrag) return;
    if (altTemplateSize > bparams.maxFrag) return;

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
        const pos_t fragOverlap(std::min((1+iparams.centerPos1-fragBeginRefPos), (fragEndRefPos-iparams.centerPos2)));
#ifdef DEBUG_MEGAPAIR
        log_os << __FUNCTION__ << ": frag begin/end/overlap: " << fragBeginRefPos << " " << fragEndRefPos << " " << fragOverlap << "\n";
#endif
        if (fragOverlap < pairOpt.minFragSupport) return;
    }

    SVFragmentEvidence& fragment(evidence.getSample(bparams.isTumor)[bamRead.qname()]);

    SVFragmentEvidenceRead& evRead(fragment.getRead(bamRead.is_first()));
    setReadEvidence(iparams.minMapQ, bamRead, evRead);

    setAlleleFrag(*bparams.fragDistroPtr, altTemplateSize, fragment.alt.getBp(isBp1));

    // when an alt entry is made for a fragment, we /*always*/ create corresponding ref entry
    // in theory this will get picked up by the ref scanner anyway, but the cost of missing this
    // is all sorts of really bad somatic FNs
    setAlleleFrag(*bparams.fragDistroPtr, templateSize, fragment.ref.getBp(isBp1));
}
