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

#include "SVCandidateAssemblyRefiner.hh"

#include "alignment/AlignmentUtil.hh"
#include "blt_util/log.hh"
#include "blt_util/samtools_fasta_util.hh"
#include "blt_util/seq_util.hh"
#include "manta/SVCandidateUtil.hh"
#include "manta/SVLocusAssembler.hh"
#include "manta/SVReferenceUtil.hh"

#include "boost/foreach.hpp"

#include <iostream>

//#define DEBUG_REFINER



/// process assembly/align info into simple reference coordinates that can be reported in the output vcf:
///
/// \param[in] isAlign1 if true, this breakend was aligned first by the jump alinger, and therefore left-aligned (if fwd) or right-aligned (if rev)
/// \param[in] jumpRange homologous range across the breakend
///
static
void
adjustAssembledBreakend(
    const Alignment& align,
    const bool isAlign1,
    const unsigned jumpRange,
    const reference_contig_segment& ref,
    const bool isReversed,
    SVBreakend& bp)
{
    const pos_t bpBeginOffset(getAlignBeginOffset(align, ref.seq().size(), isReversed));
    const pos_t bpEndOffset(getAlignEndOffset(align, ref.seq().size(), isReversed));

    const bool isBpAtAlignEnd(bp.state == SVBreakendState::RIGHT_OPEN);

    const pos_t bpBreakendOffset(isBpAtAlignEnd ? (bpEndOffset -1) : bpBeginOffset );
    const pos_t bpBreakendPos(ref.get_offset() + bpBreakendOffset);

    const bool isLeftAligned(isAlign1 == isBpAtAlignEnd);

    known_pos_range2& range(bp.interval.range);

    if (isLeftAligned)
    {
        range.set_begin_pos(bpBreakendPos);
        range.set_end_pos(bpBreakendPos + static_cast<pos_t>(jumpRange) + 1);
    }
    else
    {
        range.set_begin_pos(bpBreakendPos - static_cast<pos_t>(jumpRange));
        range.set_end_pos(bpBreakendPos + 1);
    }
}


static
bool
isFilterSpanningAlignment(
    const GlobalJumpAligner<int> aligner,
    const ALIGNPATH::path_t& apath,
    const bool isFirstRead)
{
    // require min length of each contig sub-alignment even after off-reference clipping:
    static const unsigned minAlignReadLength(30);

    // require min fraction of optimal score in each contig sub-alignmnet:
    static const float minScoreFrac(0.75);

    const unsigned readSize(apath_read_length(apath));
    const unsigned clipSize(isFirstRead ?
                            apath_soft_clip_lead_size(apath) :
                            apath_soft_clip_trail_size(apath));

    assert(clipSize <= readSize);

    const unsigned clippedReadSize(readSize-clipSize);

    if (clippedReadSize < minAlignReadLength)
    {
#ifdef DEBUG_REFINER
        log_os << "Rejecting highest scoring contig sub-alignment. isFirst?: " << isFirstRead << ". Sub-alignment read length after clipping is: " << clippedReadSize << " min size is: " << minAlignReadLength << "\n";
#endif
        return true;
    }

    int nonClipScore(aligner.getPathScore(apath, false));
    const int optimalScore(clippedReadSize * aligner.getScores().match);

    assert(optimalScore>0);
    if (nonClipScore < 0) nonClipScore = 0;

    const float scoreFrac(static_cast<float>(nonClipScore)/static_cast<float>(optimalScore));

    if (scoreFrac < minScoreFrac)
    {
#ifdef DEBUG_REFINER
        log_os << "Rejecting highest scoring contig sub-alignment. isFirst?: " << isFirstRead << ". Fraction of optimal alignment score is: " << scoreFrac << " minScoreFrac: " << minScoreFrac << "\n";
#endif
        return true;
    }
    return false;
}



/// identify path indel sequences with an insert or delete segment greater than minSize
///
static
void
getLargeIndelSegments(
    const ALIGNPATH::path_t& apath,
    const unsigned minSize,
    std::vector<std::pair<unsigned,unsigned> >& segments)
{
    using namespace ALIGNPATH;

    bool isInSegment(false);
    bool isCandidate(false);
    unsigned segmentStart(0);

    const unsigned as(apath.size());
    for (unsigned i(0); i<as; ++i)
    {
        const path_segment& ps(apath[i]);

        if ((ps.type == DELETE) || (ps.type == INSERT))
        {
            if (ps.length>=minSize) isCandidate=true;
            if (! isInSegment) segmentStart = i;
            isInSegment=true;
        }
        else
        {
            if (isCandidate)
            {
                assert(i>0);
                segments.push_back(std::make_pair(segmentStart,(i-1)));
            }
            isInSegment=false;
            isCandidate=false;
        }
    }

    if (isCandidate)
    {
        assert(as>0);
        segments.push_back(std::make_pair(segmentStart,(as-1)));
    }
}



/// add simple cigar string to spanning alignments for the subset of cases (insertions and deletions) where this is possible
///
/// note that we may not always print this out, even though we compute the cigar here -- this is dependent on the output file
/// format and conventions related to variant size, precision, etc.
///
static
void
addCigarToSpanningAlignment(
    SVCandidate& sv)
{
    const SV_TYPE::index_t svType(getSVType(sv));

    if (svType != SV_TYPE::INDEL) return;

    const bool isBp1First(sv.bp1.interval.range.begin_pos()<=sv.bp2.interval.range.begin_pos());

    const SVBreakend& bpA(isBp1First ? sv.bp1 : sv.bp2);
    const SVBreakend& bpB(isBp1First ? sv.bp2 : sv.bp1);

    const unsigned deleteSize(bpB.interval.range.begin_pos() - bpA.interval.range.begin_pos());
    const unsigned insertSize(sv.insertSeq.size());

    if(insertSize)
    {
        sv.insertAlignment.push_back(ALIGNPATH::path_segment(ALIGNPATH::INSERT,insertSize));
    }

    if(deleteSize)
    {
        sv.insertAlignment.push_back(ALIGNPATH::path_segment(ALIGNPATH::DELETE,deleteSize));
    }
}



static
bool
isSmallSVSegmentFilter(
    const AlignerBase<int>& aligner,
    const ALIGNPATH::path_t& apath,
    const bool isLeadingPath)
{
    static const unsigned minAlignRefSpan(30); ///< min reference lenght for alignment
    static const unsigned minAlignReadLength(30); ///< min length of alignment after off-reference clipping
    static const float minScoreFrac(0.75); ///< min fraction of optimal score in each contig sub-alignment:


    const unsigned refSize(apath_read_length(apath));

    if (refSize < minAlignRefSpan)
    {
        return true;
    }

    const unsigned pathSize(apath_read_length(apath));
    const unsigned clipSize(isLeadingPath ?
                            apath_soft_clip_lead_size(apath) :
                            apath_soft_clip_trail_size(apath));

    assert(clipSize <= pathSize);

    const unsigned clippedPathSize(pathSize-clipSize);

    if (clippedPathSize < minAlignReadLength)
    {
#ifdef DEBUG_REFINER
//        log_os << "Rejecting highest scoring contig sub-alignment. isFirst?: " << isFirstRead << ". Sub-alignmnet read length after clipping is: " << clippedPathSize << " min size is: " << minAlignReadLength << "\n";
#endif
        return true;
    }

    int nonClipScore(aligner.getPathScore(apath, false));
    const int optimalScore(clippedPathSize * aligner.getScores().match);

    assert(optimalScore>0);
    if (nonClipScore < 0) nonClipScore = 0;

    const float scoreFrac(static_cast<float>(nonClipScore)/static_cast<float>(optimalScore));

    if (scoreFrac < minScoreFrac)
    {
#ifdef DEBUG_REFINER
//        log_os << "Rejecting highest scoring contig sub-alignment. isFirst?: " << isFirstRead << ". Fraction of optimal alignment score is: " << scoreFrac << " minScoreFrac: " << minScoreFrac << "\n";
#endif
        return true;
    }

    return false;
}



/// test whether this single-node assembly is (1) an interesting variant above the minimum size and
/// (2) passes QC otherwise (appropriate flanking regions, etc)
///
static
bool
isFilterSmallSVAlignment(
    const GlobalAligner<int> aligner,
    const ALIGNPATH::path_t& apath,
    const unsigned minCandidateIndelSize,
    std::vector<std::pair<unsigned,unsigned> >& candidateSegments)
{
    using namespace ALIGNPATH;

    // (1) identify all indels above minimum size:
    //
    getLargeIndelSegments(apath, minCandidateIndelSize, candidateSegments);

    // escape if this is a reference or small indel alignment
    if (candidateSegments.empty()) return true;

    // test quality of alignmnet segments surrounding the variant region:
    const unsigned firstCandIndelSegment(candidateSegments.front().first);
    const unsigned lastCandIndelSegment(candidateSegments.back().second);

    const path_t leadingPath(apath.begin(), apath.begin()+firstCandIndelSegment);
    const path_t trailingPath(apath.begin()+lastCandIndelSegment+1, apath.end());

    if (isSmallSVSegmentFilter(aligner, leadingPath, true)) return true;
    if (isSmallSVSegmentFilter(aligner, trailingPath, false)) return true;

    return false;
}



/// get the range over which an alignment element can vary with equal edit distance
///
/// \param[in] refRange range of the event (ie indel) of interest in reference coordinates
/// \param[in] readRange range of the event (ie indel) of interest in read coordinates
///
/// range coordinates are zero indexed and start at the first affected positions (so are not like vcf coordinates)
/// for instance:
////  the deletion 10M1D10M would have refRange(10,11), readRange(10,10)
////  the insertion 10M1I10M would have refRange(10,10), readRange(10,11)
///
static
known_pos_range2
getVariantRange(
    const std::string& ref,
    const known_pos_range2& refRange,
    const std::string& read,
    const known_pos_range2& readRange)
{
    // check how far we can slide to the right:
    const pos_t maxRightOffset(std::min(ref.size()-refRange.end_pos(), read.size()-readRange.end_pos()));
    pos_t rightOffset(0);
    for (; rightOffset<maxRightOffset; ++rightOffset)
    {
        const char refSym(ref[refRange.begin_pos()+rightOffset]);
        const char readSym(read[readRange.begin_pos()+rightOffset]);
        if (refSym != readSym) break;
    }

    // check how far we can slide to the left:
    const pos_t minLeftOffset(std::max(-refRange.begin_pos(), -readRange.begin_pos()));
    pos_t leftOffset(0);
    for (; leftOffset>=minLeftOffset; --leftOffset)
    {
        const char refSym(ref[refRange.end_pos()+leftOffset-1]);
        const char readSym(read[readRange.end_pos()+leftOffset-1]);
        if (refSym != readSym) break;
    }

    return known_pos_range2(leftOffset,rightOffset);
}



/// process smallSV alignment section into a usable sv candidate
static
void
setSmallCandSV(
    const reference_contig_segment& ref,
    const std::string& contig,
    const Alignment& align,
    const std::pair<unsigned,unsigned>& segRange,
    SVCandidate& sv)
{
    sv.setPrecise();

    known_pos_range2 readRange;
    known_pos_range2 refRange;

    // by how many positions can the alignment position vary with the same alignment score?:
    known_pos_range2 cipos;
    {
        using namespace ALIGNPATH;

        pos_t readPos(0);
        pos_t refPos(align.beginPos);

        const path_t& apath(align.apath);
        const unsigned as(apath.size());
        for (unsigned i(0); i<as; ++i)
        {
            const path_segment& ps(apath[i]);
            if (i == segRange.first)
            {
                refRange.set_begin_pos(refPos);
                readRange.set_begin_pos(readPos);
            }

            if (is_segment_type_ref_length(ps.type)) refPos += ps.length;
            if (is_segment_type_read_length(ps.type)) readPos += ps.length;

            if (i == segRange.second)
            {
                refRange.set_end_pos(refPos);
                readRange.set_end_pos(readPos);
            }
        }

        cipos = getVariantRange(ref.seq(),refRange, contig, readRange);

        assert(cipos.begin_pos() == 0);
    }

    sv.bp1.state = SVBreakendState::RIGHT_OPEN;
    const pos_t beginPos(ref.get_offset()+refRange.begin_pos()-1);
    sv.bp1.interval.range.set_range(beginPos,beginPos+cipos.end_pos()+1);

    sv.bp2.state = SVBreakendState::LEFT_OPEN;
    const pos_t endPos(ref.get_offset()+refRange.end_pos());
    sv.bp2.interval.range.set_range(endPos,endPos+cipos.end_pos()+1);
    sv.bp2.interval.tid = sv.bp1.interval.tid;

    sv.insertSeq = contig.substr(readRange.begin_pos(),readRange.size());

    // add CIGAR for all indels:
    sv.insertAlignment = ALIGNPATH::path_t(align.apath.begin()+segRange.first, align.apath.begin()+segRange.second+1);
}



SVCandidateAssemblyRefiner::
SVCandidateAssemblyRefiner(
    const GSCOptions& opt,
    const bam_header_info& header) :
    _opt(opt),
    _header(header),
    _smallSVAssembler(opt.scanOpt, opt.refineOpt.smallSVAssembleOpt, opt.alignFileOpt, opt.statsFilename),
    _spanningAssembler(opt.scanOpt, opt.refineOpt.spanningAssembleOpt, opt.alignFileOpt, opt.statsFilename),
    _smallSVAligner(opt.refineOpt.smallSVAlignScores),
    _spanningAligner(opt.refineOpt.spanningAlignScores, opt.refineOpt.jumpScore)
{}



void
SVCandidateAssemblyRefiner::
getCandidateAssemblyData(
    const SVCandidate& sv,
    const SVCandidateSetData& /*svData*/,
    SVCandidateAssemblyData& assemblyData) const
{
#ifdef DEBUG_REFINER
    static const std::string logtag("getCandidateAssemblyData");
    log_os << logtag << " START sv: " << sv;
#endif

    assemblyData.clear();

    // separate the problem into different assembly categories:
    //
    if (isSimpleBreakend(sv.bp1.state) && isSimpleBreakend(sv.bp1.state))
    {
        // this case assumes two suspected breakends with a direction to each, most common large scale SV case:
        getJumpAssembly(sv, assemblyData);
    }
    else if ((sv.bp1.state == SVBreakendState::COMPLEX))
    {
        // this case assumes a single-interval local assembly, this is the most common case for small-scale SVs/indels
        getSmallSVAssembly(sv, assemblyData);
    }
    else
    {
        log_os << "Unknown candidate SV: " << sv << "\n";
        assert(false && "Unknown candidate SV type");
    }
}



void
SVCandidateAssemblyRefiner::
getJumpAssembly(
    const SVCandidate& sv,
    SVCandidateAssemblyData& assemblyData) const
{
#ifdef DEBUG_REFINER
    static const std::string logtag("getJumpAssembly");
    log_os << logtag << " START\n";
#endif

    // how much additional reference sequence should we extract from around
    // each side of the breakend region?
    static const pos_t extraRefEdgeSize(200);

    // if the breakends have a simple insert/delete orientation and the alignment regions overlap, then handle this case as
    // a local assembly problem:
    if (sv.bp1.interval.tid == sv.bp2.interval.tid)
    {
        if (! SVBreakendState::isSameOrientation(sv.bp1.state,sv.bp2.state))
        {
            const SV_TYPE::index_t svType(getSVType(sv));
            if ((svType == SV_TYPE::INDEL) || (svType == SV_TYPE::COMPLEX))
            {
                if ( isRefRegionOverlap( _header, extraRefEdgeSize, sv) )
                {
                    // transform SV into a single region format:
                    SVCandidate singleSV = sv;
                    singleSV.bp1.state = SVBreakendState::COMPLEX;
                    singleSV.bp2.state = SVBreakendState::UNKNOWN;
                    singleSV.bp1.interval.range.merge_range(sv.bp2.interval.range);

#ifdef DEBUG_REFINER
                    log_os << logtag << " candidate breakends regions are too close, transferring problem to local assembler\n";
#endif

                    getSmallSVAssembly(singleSV, assemblyData);
                    return;
                }
            }
        }
    }

    assemblyData.isSpanning = true;

    //
    // based on sv candidate, we classify the expected relationship between the contig and the sv breakends:
    //
    bool isBp2AlignedFirst(false); ///< should the contig on the fwd strand align bp2->bp1 (true) or bp1->bp2 (false)

    bool isBp1Reversed(false); ///< should all bp1 reads be reversed for the contig to assemble correctly?
    bool isBp2Reversed(false); ///< should all bp2 reads be reversed for the contig to assemble correctly?

    if (sv.bp1.state != sv.bp2.state)
    {
        // if there's one right-open breakend and one left-open breakend, no matter the bp1/bp2 chromosome and
        // relative bp1/bp2 order etc. we:
        // 1. don't need to do any read/reference reversals
        // 2. always treat the right-open breakend as the first alignment region in order:
        //
        if (sv.bp2.state == SVBreakendState::RIGHT_OPEN)
        {
            isBp2AlignedFirst = true;
        }
    }
    else
    {
        // If both breakends open in the same direction, then:
        // 1. the reads from one breakend need to be reversed
        // 2. the reference from that same breakend needs to be reversed
        // 3. Treat the un-reversed RIGHT_OPEN or reversed LEFT_OPEN as the first alignment region in order
        //      Note that in the scheme below, we chose which bp to reverse so that no-reordering is required
        //
        if (sv.bp1.state == SVBreakendState::RIGHT_OPEN)
        {
            isBp2Reversed = true;
        }
        else
        {
            isBp1Reversed = true;
        }
    }

    assemblyData.isBp2AlignedFirst = isBp2AlignedFirst;
    assemblyData.isBp1Reversed = isBp1Reversed;
    assemblyData.isBp2Reversed = isBp2Reversed;

    // assemble contig spanning the breakend:
    _spanningAssembler.assembleSVBreakends(sv.bp1, sv.bp2, isBp1Reversed, isBp2Reversed, assemblyData.contigs);

    getSVReferenceSegments(_opt.referenceFilename, _header, extraRefEdgeSize, sv, assemblyData.bp1ref, assemblyData.bp2ref);
    std::string bp1refSeq = assemblyData.bp1ref.seq();
    std::string bp2refSeq = assemblyData.bp2ref.seq();
    if (isBp1Reversed) reverseCompStr(bp1refSeq);
    if (isBp2Reversed) reverseCompStr(bp2refSeq);
    const std::string* align1RefStrPtr(&bp1refSeq);
    const std::string* align2RefStrPtr(&bp2refSeq);

    if (isBp2AlignedFirst) std::swap(align1RefStrPtr, align2RefStrPtr);

#ifdef DEBUG_REFINER
    log_os << logtag << " al1Ref: " << *align1RefStrPtr << "\n";
    log_os << logtag << " al2Ref: " << *align2RefStrPtr << "\n";
#endif

    const unsigned contigCount(assemblyData.contigs.size());

#ifdef DEBUG_REFINER
    log_os << logtag << " contigCount: " << contigCount << "\n";
    for (unsigned contigIndex(0); contigIndex<contigCount; ++contigIndex)
    {
        const AssembledContig& contig(assemblyData.contigs[contigIndex]);
        log_os << logtag << " contigIndex: " << contigIndex << " contig: " << contig;
    }
#endif

    // make sure an alignment object exists for every contig, even if it's empty:
    assemblyData.spanningAlignments.resize(contigCount);

    bool isHighScore(false);
    unsigned highScoreIndex(0);

    for (unsigned contigIndex(0); contigIndex<contigCount; ++contigIndex)
    {
        const AssembledContig& contig(assemblyData.contigs[contigIndex]);

#ifdef DEBUG_REFINER
        log_os << logtag << " start aligning contigIndex: " << contigIndex << "\n";
#endif

        JumpAlignmentResult<int>& alignment(assemblyData.spanningAlignments[contigIndex]);

        _spanningAligner.align(
            contig.seq.begin(), contig.seq.end(),
            align1RefStrPtr->begin(), align1RefStrPtr->end(),
            align2RefStrPtr->begin(), align2RefStrPtr->end(),
            alignment);

        std::string extendedContig;
        getExtendedContig(alignment, contig.seq, align1RefStrPtr, align2RefStrPtr, extendedContig);
        assemblyData.extendedContigs.push_back(extendedContig);

#ifdef DEBUG_REFINER
        log_os << logtag << " contigIndex: " << contigIndex << " alignment: " << alignment;

        std::string bp1Seq,bp2Seq,insertSeq;
        getFwdStrandQuerySegments(alignment, contig.seq,
                                  isBp2AlignedFirst, isBp1Reversed, isBp2Reversed,
                                  bp1Seq, bp2Seq, insertSeq);
        log_os << logtag << "\tbp1seq_fwd: " << bp1Seq << "\n";
        log_os << logtag << "\tinsseq_fwd: " << insertSeq << "\n";
        log_os << logtag << "\tbp2seq_fwd: " << bp2Seq << "\n";
#endif

        // QC the alignment to make sure it spans the two breakend locations:
        static const unsigned minAlignRefSpan(20);
        const bool isAlignment1Good(alignment.align1.isAligned() && (apath_ref_length(alignment.align1.apath) >= minAlignRefSpan));
        const bool isAlignment2Good(alignment.align2.isAligned() && (apath_ref_length(alignment.align2.apath) >= minAlignRefSpan));
        const bool isAlignmentGood(isAlignment1Good && isAlignment2Good);

        if (! isAlignmentGood) continue;

        if ((! isHighScore) || (alignment.score > assemblyData.spanningAlignments[highScoreIndex].score))
        {
            isHighScore = true;
            highScoreIndex=contigIndex;
        }
    }

    if (! isHighScore) return;

    // set any additional QC steps before deciding an alignment is usable:

    // check the min size and fraction of optimal score for each sub-alignment:
    {
        const SVCandidateAssemblyData::JumpAlignmentResultType& hsAlign(assemblyData.spanningAlignments[highScoreIndex]);

        if (isFilterSpanningAlignment(_spanningAligner, hsAlign.align1.apath, true)) return;
        if (isFilterSpanningAlignment(_spanningAligner, hsAlign.align2.apath, false)) return;
    }

    // TODO: min context, etc.


    // ok, passed QC -- mark the high-scoring alignment as usable for hypothesis refinement:
    {
        assemblyData.bestAlignmentIndex = highScoreIndex;
#ifdef DEBUG_REFINER
        log_os << logtag << " highscoreid: " << highScoreIndex << " alignment: " << assemblyData.spanningAlignments[highScoreIndex];
#endif

        // process the alignment into information that's easily usable in the vcf output
        // (ie. breakends in reference coordinates)

        const AssembledContig& bestContig(assemblyData.contigs[assemblyData.bestAlignmentIndex]);
        const SVCandidateAssemblyData::JumpAlignmentResultType& bestAlign(assemblyData.spanningAlignments[assemblyData.bestAlignmentIndex]);

        // first get each alignment associated with the correct breakend:
        const Alignment* bp1AlignPtr(&bestAlign.align1);
        const Alignment* bp2AlignPtr(&bestAlign.align2);

        if (isBp2AlignedFirst) std::swap(bp1AlignPtr, bp2AlignPtr);

        // summarize usable output information in a second SVBreakend object -- this is the 'refined' sv:
        assemblyData.svs.push_back(sv);
        SVCandidate& newSV(assemblyData.svs.back());
        newSV.assemblyIndex = 0;

        newSV.setPrecise();

        adjustAssembledBreakend(*bp1AlignPtr, (! isBp2AlignedFirst), bestAlign.jumpRange, assemblyData.bp1ref, isBp1Reversed, newSV.bp1);
        adjustAssembledBreakend(*bp2AlignPtr, (isBp2AlignedFirst), bestAlign.jumpRange, assemblyData.bp2ref, isBp2Reversed, newSV.bp2);

        // fill in insertSeq:
        newSV.insertSeq.clear();
        if (bestAlign.jumpInsertSize > 0)
        {
            getFwdStrandInsertSegment(bestAlign, bestContig.seq, isBp1Reversed, newSV.insertSeq);
        }

        // add CIGAR for any simple (insert/delete) cases:
        addCigarToSpanningAlignment(newSV);

#ifdef DEBUG_REFINER
        log_os << logtag << " highscore refined sv: " << newSV;
#endif
    }
}



void
SVCandidateAssemblyRefiner::
getSmallSVAssembly(
    const SVCandidate& sv,
    SVCandidateAssemblyData& assemblyData) const
{
#ifdef DEBUG_REFINER
    static const std::string logtag("getSmallSVAssembly");
    log_os << logtag << " START\n";
#endif

    assemblyData.isSpanning = false;

    // assemble contigs in the breakend region
    _smallSVAssembler.assembleSingleSVBreakend(sv.bp1, assemblyData.contigs);

    // how much additional reference sequence should we extract from around
    // each side of the breakend region?
    static const pos_t extraRefEdgeSize(700);

    // min alignment context
    //const unsigned minAlignContext(4);

    getIntervalReferenceSegment(_opt.referenceFilename, _header, extraRefEdgeSize, sv.bp1.interval, assemblyData.bp1ref);
    const std::string* align1RefStrPtr(&assemblyData.bp1ref.seq());

#ifdef DEBUG_REFINER
    log_os << logtag << " al1Ref: " << *align1RefStrPtr << "\n";
#endif

    const unsigned contigCount(assemblyData.contigs.size());

#ifdef DEBUG_REFINER
    log_os << logtag << " contigCount: " << contigCount << "\n";
    for (unsigned contigIndex(0); contigIndex<contigCount; ++contigIndex)
    {
        const AssembledContig& contig(assemblyData.contigs[contigIndex]);
        log_os << logtag << " contigIndex: " << contigIndex << " contig: " << contig;
    }
#endif

    // make sure an alignment object exists for every contig, even if it's empty:
    assemblyData.smallSVAlignments.resize(contigCount);
    assemblyData.smallSVSegments.resize(contigCount);

    bool isHighScore(false);
    unsigned highScoreIndex(0);

    for (unsigned contigIndex(0); contigIndex<contigCount; ++contigIndex)
    {
        const AssembledContig& contig(assemblyData.contigs[contigIndex]);

#ifdef DEBUG_REFINER
        log_os << logtag << " start aligning contigIndex: " << contigIndex << "\n";
#endif

        AlignmentResult<int>& alignment(assemblyData.smallSVAlignments[contigIndex]);

        _smallSVAligner.align(
            contig.seq.begin(), contig.seq.end(),
            align1RefStrPtr->begin(), align1RefStrPtr->end(),
            alignment);

        // remove candidate from consideration unless we find a sufficiently large indel with good flanking sequence:
        std::vector<std::pair<unsigned,unsigned> >& candidateSegments(assemblyData.smallSVSegments[contigIndex]);
        const bool isFilterSmallSV( isFilterSmallSVAlignment(_smallSVAligner, alignment.align.apath, _opt.scanOpt.minCandidateVariantSize, candidateSegments));

#ifdef DEBUG_REFINER
        log_os << logtag << " contigIndex: " << contigIndex << " isFilter " << isFilterSmallSV << " alignment: " << alignment;
#endif

        if (isFilterSmallSV) continue;

        // keep the highest scoring QC'd candidate:
        // TODO: we should keep all QC'd candidates for the small event case
        // FIXME : prevents us from finding overlapping events, keep vector of high-scoring contigs?
        if ((! isHighScore) || (alignment.score > assemblyData.smallSVAlignments[highScoreIndex].score))
        {
#ifdef DEBUG_REFINER
            log_os << logtag << " contigIndex: " << contigIndex << " is high score\n";
#endif
            isHighScore = true;
            highScoreIndex=contigIndex;
        }
    }

    if (! isHighScore) return;

    // set any additional QC steps before deciding an alignment is usable:
    // TODO:


    // ok, passed QC -- mark the high-scoring alignment as usable for hypothesis refinement:
    {
        assemblyData.bestAlignmentIndex = highScoreIndex;
#ifdef DEBUG_REFINER
        log_os << logtag << " highscoreid: " << highScoreIndex << " alignment: " << assemblyData.smallSVAlignments[highScoreIndex];
#endif

        // process the alignment into information that's easily usable in the vcf output
        // (ie. breakends in reference coordinates)

        const AssembledContig& bestContig(assemblyData.contigs[assemblyData.bestAlignmentIndex]);
        const SVCandidateAssemblyData::SmallAlignmentResultType& bestAlign(assemblyData.smallSVAlignments[assemblyData.bestAlignmentIndex]);

        const SVCandidateAssemblyData::CandidateSegmentSetType& candidateSegments(assemblyData.smallSVSegments[assemblyData.bestAlignmentIndex]);
        BOOST_FOREACH(const SVCandidateAssemblyData::CandidateSegmentType& segRange, candidateSegments)
        {
            assemblyData.svs.push_back(sv);
            SVCandidate& newSV(assemblyData.svs.back());
            newSV.assemblyIndex = (assemblyData.svs.size() - 1);
            setSmallCandSV(assemblyData.bp1ref, bestContig.seq, bestAlign.align, segRange, newSV);

#ifdef DEBUG_REFINER
            log_os << logtag << "small refined sv: " << newSV;
#endif
        }
    }
}
