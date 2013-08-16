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

#include "SVCandidateAssemblyRefiner.hh"

#include "alignment/AlignmentUtil.hh"
#include "blt_util/samtools_fasta_util.hh"
#include "blt_util/seq_util.hh"
#include "manta/SVLocusAssembler.hh"
#include "manta/SVReferenceUtil.hh"

#include "boost/foreach.hpp"

#define DEBUG_REFINER

#ifdef DEBUG_REFINER
#include <iostream>
#include "blt_util/log.hh"
#endif


/// process assembly/align info into simple reference coordinates that can be reported in the output vcf:
static
void
adjustAssembledBreakend(
    const Alignment& align,
    const reference_contig_segment& ref,
    const bool isReversed,
    SVBreakend& bp)
{
    const pos_t bpBeginOffset(getAlignBeginOffset(align, ref.seq().size(), isReversed));
    const pos_t bpEndOffset(getAlignEndOffset(align, ref.seq().size(), isReversed));

    const bool isBpAtAlignEnd(bp.state == SVBreakendState::RIGHT_OPEN);

    const pos_t bpBreakendOffset(isBpAtAlignEnd ? bpEndOffset : bpBeginOffset+1);
    const pos_t bpBreakendPos(ref.get_offset() + bpBreakendOffset);

    known_pos_range2& range(bp.interval.range);
    range.set_begin_pos(bpBreakendPos);
    range.set_end_pos(bpBreakendPos);
}



SVCandidateAssemblyRefiner::
SVCandidateAssemblyRefiner(
    const GSCOptions& opt,
    const bam_header_info& header,
    const SmallAssemblerOptions& assembleOpt,
    const AlignmentScores<int>& alignScore,
    const int jumpScore) :
    _opt(opt),
    _svAssembler(opt, assembleOpt),
    _header(header),
    _aligner(alignScore,jumpScore)
{}



void
SVCandidateAssemblyRefiner::
getCandidateAssemblyData(
    const SVCandidate& sv,
    const SVCandidateSetData& /*svData*/,
    SVCandidateAssemblyData& adata) const
{
#ifdef DEBUG_REFINER
    log_os << "getCandidateAssemblyData START sv: " << sv;
#endif

    adata.clear();

    // filter out simple single region breakends for now:
    if (! (isSimpleBreakend(sv.bp1.state) && isSimpleBreakend(sv.bp1.state))) return;

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

    // assemble contig spanning the breakend:
    _svAssembler.assembleSVBreakends(sv.bp1, sv.bp2, isBp1Reversed, isBp2Reversed, adata.contigs);


    // how much additional reference sequence should we extract from around
    // each side of the breakend region?
    static const pos_t extraRefEdgeSize(300);

    // min alignment context
    //const unsigned minAlignContext(4);
    // don't align contigs shorter than this
    static const unsigned minContigLen(75);

    reference_contig_segment bp1ref, bp2ref;
    getSVReferenceSegments(_opt.referenceFilename, _header, extraRefEdgeSize, sv, bp1ref, bp2ref);
    const std::string* align1RefStrPtr(&bp1ref.seq());
    const std::string* align2RefStrPtr(&bp2ref.seq());

    if (isBp1Reversed) reverseCompStr(bp1ref.seq());
    if (isBp2Reversed) reverseCompStr(bp2ref.seq());

    if (isBp2AlignedFirst) std::swap(align1RefStrPtr, align2RefStrPtr);

    const unsigned contigCount(adata.contigs.size());

#ifdef DEBUG_REFINER
    log_os << "contigCount: " << contigCount << "\n";
    for (unsigned contigIndex(0); contigIndex<contigCount; ++contigIndex)
    {
        const AssembledContig& contig(adata.contigs[contigIndex]);
        log_os << "cid: " << contigIndex << " contig: " << contig;
    }
#endif

    // make sure an alignment object exists for every contig, even if it's empty:
    adata.alignments.resize(contigCount);

    bool isHighScore(false);
    unsigned highScoreIndex(0);

    for (unsigned contigIndex(0); contigIndex<contigCount; ++contigIndex)
    {
        const AssembledContig& contig(adata.contigs[contigIndex]);

        // QC contig prior to alignment:
        if (contig.seq.size() < minContigLen) continue;

        JumpAlignmentResult<int>& alignment(adata.alignments[contigIndex]);

        _aligner.align(
            contig.seq.begin(), contig.seq.end(),
            align1RefStrPtr->begin(), align1RefStrPtr->end(),
            align2RefStrPtr->begin(), align2RefStrPtr->end(),
            alignment);

#ifdef DEBUG_REFINER
        log_os << "cid: " << contigIndex << " alignment: " << alignment;
#endif

        // QC the alignment to make sure it spans the two breakend locations:
        static const unsigned minAlignRefSpan(20);
        const bool isAlignment1Good(alignment.align1.isAligned() && (apath_ref_length(alignment.align1.apath) >= minAlignRefSpan));
        const bool isAlignment2Good(alignment.align2.isAligned() && (apath_ref_length(alignment.align2.apath) >= minAlignRefSpan));
        const bool isAlignmentGood(isAlignment1Good && isAlignment2Good);

        if (! isAlignmentGood) continue;

        if ((! isHighScore) || (alignment.score > adata.alignments[highScoreIndex].score))
        {
            isHighScore = true;
            highScoreIndex=contigIndex;
        }
    }

    // set any additional QC steps before deciding an alignment is usable:

    // TODO: min context, etc.


    // ok, passed QC -- mark the high-scoring alignment as usable for hypothesis refinement:
    if (isHighScore)
    {
        adata.isBestAlignment = true;
        adata.bestAlignmentIndex = highScoreIndex;
#ifdef DEBUG_REFINER
        log_os << "highscoreid: " << highScoreIndex << " alignment: " << adata.alignments[highScoreIndex];
#endif

        // process the alignment into information that's easily usable in the vcf output
        // (ie. breakends in reference coordinates)

        // copy the best alignment so that we can revese coordinates,etc:
        SVCandidateAssemblyData::JumpAlignmentResultType bestAlignCopy(adata.alignments[adata.bestAlignmentIndex]);

        // first get each alignment associated with the correct breakend:
        Alignment* bp1AlignPtr(&bestAlignCopy.align1);
        Alignment* bp2AlignPtr(&bestAlignCopy.align2);

        if (isBp2AlignedFirst) std::swap(bp1AlignPtr, bp2AlignPtr);

        adata.sv = sv;

        adjustAssembledBreakend(*bp1AlignPtr,bp1ref,isBp1Reversed,adata.sv.bp1);
        adjustAssembledBreakend(*bp2AlignPtr,bp2ref,isBp2Reversed,adata.sv.bp2);

#ifdef DEBUG_REFINER
        log_os << "highscore refined sv: " << adata.sv;
#endif
    }
}
