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
#include "blt_util/seq_util.hh"
#include "blt_util/SimpleAlignment.hh"
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



ContigParams::
ContigParams(
    const SVCandidateAssemblyData& assemblyData,
    const SVCandidate& sv) :
    extSeq(assemblyData.extendedContigs[sv.assemblyAlignIndex]),
    beginPos(assemblyData.bp1ref.get_offset()),
    endPos(beginPos + assemblyData.bp1ref.seq().size())
{
    // get offsets of breakpoints in the extended contig
    const AlignmentResult<int>& alignment(assemblyData.smallSVAlignments[sv.assemblyAlignIndex]);
    const std::pair<unsigned, unsigned>& alignSegment(assemblyData.smallSVSegments[sv.assemblyAlignIndex][sv.assemblySegmentIndex]);

    const pos_t bp1HomLength(sv.bp1.interval.range.size()-1);
    const pos_t bp2HomLength(sv.bp2.interval.range.size()-1);
    assert(bp1HomLength >= 0);
    assert(bp2HomLength >= 0);

    ALIGNPATH::path_t apathTillSvStart(&alignment.align.apath[0], &alignment.align.apath[alignSegment.first]);
    ALIGNPATH::path_t apathTillSvEnd(&alignment.align.apath[0], &alignment.align.apath[alignSegment.second+1]);

    // the beginPos of align is the length of reference padding in the extended contig
    // |ref padding| + |alignment segments|
    // both bp1 and bp2 include the insert and homology,
    // which can avoid false split-read evidence from normal sample when the homology is long

    /// all offset range 'begin' values correspond to the zero-indexed base immediately before the breakend on the fwd-strand,
    /// and 'end' values correspond to the zero-indexed base immediately before the breakend on the forward strand+homology range
    /// In the absence of homology, begin and end should be equal.

    /// note that we add align.beginPos here to reflect coordinates in the extended Contig, in the regular contig we wouldn't add this
    bp1Offset.set_begin_pos(alignment.align.beginPos + apath_read_length(apathTillSvStart) - 1);
    bp1Offset.set_end_pos(bp1Offset.begin_pos() + bp1HomLength);
    bp2Offset.set_begin_pos(alignment.align.beginPos + apath_read_length(apathTillSvEnd) - 1);
    bp2Offset.set_end_pos(bp2Offset.begin_pos() + bp2HomLength);
}



void
SVScorePairAltProcessor::
checkInput(
    const SVCandidateAssemblyData& assemblyData,
    const SVCandidate& sv)
{
    using namespace illumina::common;

    // this test is a bit redundant with those below, but an ounce of prevention...
    const bool isSpanning(assemblyData.isSpanning);
    assert(! isSpanning);

    // this class is designed for simple alts only:
    assert(sv.bp1.interval.tid == sv.bp2.interval.tid);
    assert(getSVType(sv) == SV_TYPE::INDEL);

    /// In case of breakend microhomology approximate the breakend as a point event at the center of the possible range:
    const pos_t centerPos1(sv.bp1.interval.range.center_pos());
    const pos_t centerPos2(sv.bp2.interval.range.center_pos());
    if (centerPos2 <= centerPos1)
    {
        std::ostringstream oss;
        oss << "ERROR: Unexpected breakend orientation in pair support routine for sv: " << sv << "\n";
        BOOST_THROW_EXCEPTION(LogicException(oss.str()));
    }
}



/// test whether a frag reference span provides sufficient support for a breakpoint of this sv:
bool
SVScorePairAltProcessor::
testFragOverlap(
    const int fragBeginRefPos,
    const int fragEndRefPos) const
{
    const pos_t fragOverlap(std::min((1+svParams.centerPos1-fragBeginRefPos), (fragEndRefPos-svParams.centerPos2)));
#ifdef DEBUG_MEGAPAIR
    log_os << __FUNCTION__ << ": frag begin/end/overlap: " << fragBeginRefPos << " " << fragEndRefPos << " " << fragOverlap << "\n";
#endif
    return (fragOverlap >= pairOpt.minFragSupport);
}



bool
SVScorePairAltProcessor::
alignShadowRead(
    const bam_record& bamRead,
    int& altTemplateSize)
{
    // TODO: basecall qualities??

    // does the shadow occur to the left or right of the insertion?
    const bool isLeftOfInsert(bamRead.is_mate_fwd_strand());

    // do we need to revcomp the sequence?
    std::string read(bamRead.get_bam_read().get_string());
    if (isLeftOfInsert)
    {
        reverseCompStr(read);
    }

    // TODO: break the contig into two parts for incomplete insertions
    AlignmentResult<int> readAlignment;

    typedef std::string::const_iterator siter;
    const siter readBegin(read.begin());
    const siter readEnd(read.end());
    siter contigBegin(_contig.extSeq.begin());
    siter contigEnd(_contig.extSeq.end());

    // if the insertion is not fully assembled, align to only part of the contig:
    int contigBeginOffset(0);
    if (sv.isUnknownSizeInsertion)
    {
        // TODO check these results in test case:
        if (isLeftOfInsert)
        {
            contigEnd = contigBegin + _contig.bp1Offset.begin_pos() + sv.unknownSizeInsertionLeftSeq.size();
        }
        else
        {
            contigBeginOffset = static_cast<int>(_contig.bp2Offset.begin_pos()) - sv.unknownSizeInsertionRightSeq.size();
            assert(contigBeginOffset>=0);
            contigBegin = contigBegin + contigBeginOffset;
        }
    }

    _shadowAligner.align(
        readBegin, readEnd,
        contigBegin, contigEnd,
        readAlignment);

    //
    // first determine if the read meets some minimal quality criteria
    //

    // require the complete alignment score to be some percentage of optimal after trimming off any expected softclip
    {
        using namespace ALIGNPATH;

        const path_t readPath(readAlignment.align.apath);

        const unsigned readSize(read.size());
        unsigned clipSize(0);

        if (sv.isUnknownSizeInsertion)
        {
            if (isLeftOfInsert)
            {
                clipSize=apath_soft_clip_trail_size(readPath);
            }
            else
            {
                clipSize=apath_soft_clip_lead_size(readPath);
            }
        }

        assert(clipSize <= readSize);

        const unsigned clippedReadSize(readSize-clipSize);

        static const unsigned minAlignReadLength(40);
        if (clippedReadSize < minAlignReadLength)
        {
            return false;
        }

        int nonClipScore(_shadowAligner.getPathScore(readPath, false));

        static const float minScoreFrac(0.85);
        const int optimalScore(clippedReadSize*_shadowAligner.getScores().match);

        const float scoreFrac(static_cast<float>(nonClipScore)/static_cast<float>(optimalScore));

        if (scoreFrac < minScoreFrac)
        {
            return false;
        }
    }

    //
    // next determine what the altTemplateSize is if we believe the alignment
    //

    known_pos_range2 fakeRefSpan;
    if (isLeftOfInsert)
    {
        fakeRefSpan.set_begin_pos(bamRead.mate_pos()-1);

        /// set fake end as if the insert allele continued in reference coordinates
        fakeRefSpan.set_end_pos(_contig.beginPos + contigBeginOffset + readAlignment.align.beginPos);
    }
    else
    {
        using namespace ALIGNPATH;
        SimpleAlignment bamAlign(bamRead);
        fakeRefSpan.set_end_pos(bamRead.mate_pos()-1 + apath_ref_length(bamAlign.path));

        /// set fake begin as if the insert allele continued TO THE LEFT in reference coordinates:
        const int readContigBeginOffset(contigBeginOffset + static_cast<int>(readAlignment.align.beginPos));
        const int readContigEndOffset(static_cast<int>(_contig.extSeq.size()) - readContigBeginOffset);
        assert(readContigEndOffset >= 0);

        fakeRefSpan.set_begin_pos(_contig.endPos - readContigEndOffset);
    }

    if (fakeRefSpan.begin_pos() > fakeRefSpan.end_pos())
    {
        using namespace illumina::common;

        std::ostringstream oss;
        oss << "ERROR: Failed to parse fragment range from bam record. Frag begin,end: " << fakeRefSpan << " bamRecord: " << bamRead << "\n";
        BOOST_THROW_EXCEPTION(LogicException(oss.str()));
    }

    altTemplateSize=fakeRefSpan.size();

    //
    // finally determine if we cross a breakend boundary
    //
    if (! testFragOverlap(fakeRefSpan.begin_pos(), fakeRefSpan.end_pos())) return false;

    // made it!

    return true;
}



void
SVScorePairAltProcessor::
processClearedRecord(
    const bam_record& bamRead)
{
    using namespace illumina::common;

    assert(bamParams.isSet);

    const pos_t refPos(bamRead.pos()-1);
    if (! bamParams.interval.range.is_pos_intersect(refPos)) return;

    // many special rules applied for large insertions:
    const bool isLargeInsert(isLargeInsertSV(sv));

    bool isShadowAlignment(false);

    int templateSize(0);
    int altTemplateSize(0);

    if (isLargeInsert)
    {
        // test for shadow
        const bool isShadowRead(_shadow.check(bamRead));

        if (isShadowRead)
        {
            // no need to check if we've encountered the shadow before b/c of left/right orientation restriction

            // don't forget to revcomp!
            /// what information do we need back form this function?
            /// 1) bool -- did the read align within the bounds of the insert region and meet some alignment quality threshold
            /// 2) if (1) provide the template size estimate:
            isShadowAlignment=alignShadowRead(bamRead,altTemplateSize);

            if (! isShadowAlignment) return;
        }
        else
        {
            // ok, not a shadow read, kick it out:
            if (bamRead.is_unmapped() || bamRead.is_mate_unmapped()) return;
        }

        // test for MAPQ0 pair
        const bool isRepeatChimera(false);
        if (isRepeatChimera)
        {

        }
        else
        {
            /// if we establish it's not a repeat chimera, then kick the candidate out:
            if (! (bamRead.is_unmapped() || bamRead.is_mate_unmapped()))
            {
                if (! is_innie_pair(bamRead)) return;
            }
        }
    }

    const bool isRealignedTemplate(isLargeInsert && isShadowAlignment);


    // for now, we allow alt evidence that aligns within the insert only, any
    // fragments that align to the reference on both sides and are simply "short"
    // due to the insertion are not considered
    /// TODO: re-evaluate this filter after shadow scoring is mature
    if (isLargeInsert && (! isShadowAlignment)) return;

#ifdef DEBUG_MEGAPAIR
    log_os << __FUNCTION__ << ": read: " << bamRead << "\n";
#endif

    /// check if fragment is too big or too small:
    if (! isRealignedTemplate)
    {
        templateSize=(std::abs(bamRead.template_size()));
        altTemplateSize=(templateSize-svParams.altShift);
    }
    if (altTemplateSize < bamParams.minFrag) return;
    if (altTemplateSize > bamParams.maxFrag) return;

    // get fragment range and check overlap with breakend:
    if (! isRealignedTemplate)
    {
        // count only from the down stream reads
        const bool isFirstBamRead(isFirstRead(bamRead));

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

        if (! testFragOverlap(fragBeginRefPos, fragEndRefPos)) return;
    }

    SVFragmentEvidence& fragment(evidence.getSample(bamParams.isTumor)[bamRead.qname()]);

    SVFragmentEvidenceRead& evRead(fragment.getRead(bamRead.is_first()));

    setReadEvidence(svParams.minMapQ, svParams.minTier2MapQ, bamRead, isShadowAlignment, evRead);

    setAlleleFrag(*bamParams.fragDistroPtr, altTemplateSize, fragment.alt.getBp(isBp1));

    if (! isRealignedTemplate)
    {
        // when an alt entry is made for a fragment, we try to always create corresponding ref entry
        // in theory this will get picked up by the ref scanner anyway, but the cost of missing this
        // is all sorts of really bad somatic FNs
        setAlleleFrag(*bamParams.fragDistroPtr, templateSize, fragment.ref.getBp(isBp1));
    }
}
