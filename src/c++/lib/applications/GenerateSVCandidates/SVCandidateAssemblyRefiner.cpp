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
    if(! (isSimpleBreakend(sv.bp1.state) && isSimpleBreakend(sv.bp1.state))) return;

    //
    // based on sv candidate, we classify the expected relationship between the contig and the sv breakends:
    //
    bool isBp2AlignedFirst(false);

    bool isBp1Reversed(false);
    bool isBp2Reversed(false);

    if(sv.bp1.state != sv.bp2.state)
    {
        // if there's one right-open breakend and one left-open breakend, no matter the chromsome, order etc. we:
        // 1. don't need to do any reversals
        // 2. always treat the right open breakend as the first alignment region in order:
        //
        if(sv.bp2.state == SVBreakendState::RIGHT_OPEN)
        {
            isBp2AlignedFirst = true;
        }
    }
    else
    {
        // if both breakends open in the same direction, then:
        // 1. the reads from one breakend need to be reversed
        // 2. the reference from that same breakend needs to be reversed
        // 3. Treat the unreversed RIGHT_OPEN or reversed LEFT_OPEN as the first alignment region in order
        //
        if(sv.bp1.state == SVBreakendState::RIGHT_OPEN)
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

    reference_contig_segment bp1ref,bp2ref;
    getSVReferenceSegments(_opt.referenceFilename, _header, extraRefEdgeSize, sv, bp1ref, bp2ref);
    const std::string* align1RefStrPtr(&bp1ref.seq());
    const std::string* align2RefStrPtr(&bp2ref.seq());

    if(isBp1Reversed) reverseCompStr(bp1ref.seq());
    if(isBp2Reversed) reverseCompStr(bp2ref.seq());

    if(isBp2AlignedFirst) std::swap(align1RefStrPtr, align2RefStrPtr);

    const unsigned contigCount(adata.contigs.size());

#ifdef DEBUG_REFINER
    log_os << "contigCount: " << contigCount << "\n";
    for(unsigned contigIndex(0);contigIndex<contigCount; ++contigIndex)
    {
        const AssembledContig& contig(adata.contigs[contigIndex]);
        log_os << "cid: " << contigIndex << " contig: " << contig;
    }
#endif

    // make sure an alignment object exists for every contig, even if it's empty:
    adata.alignments.resize(contigCount);

    bool isHighScore(false);
    unsigned highScoreIndex(0);

    for(unsigned contigIndex(0);contigIndex<contigCount; ++contigIndex)
    {
        const AssembledContig& contig(adata.contigs[contigIndex]);

        // QC contig prior to alignment:
        if (contig.seq.size() < minContigLen) continue;

        JumpAlignmentResult<int> alignment(adata.alignments[contigIndex]);

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

        if(! isAlignmentGood) continue;

        if((! isHighScore) || (alignment.score > adata.alignments[highScoreIndex].score))
        {
            isHighScore = true;
            highScoreIndex=contigIndex;
        }
    }

    // set any additional QC steps before deciding an alignment is usable:

    // TODO: min context, etc.


    // ok, passed QC -- mark the high-scoring alignment as usable:
    if(isHighScore)
    {
        adata.isBestAlignment = true;
        adata.bestAlignmentIndex = highScoreIndex;
#ifdef DEBUG_REFINER
        log_os << "highscoreid: " << highScoreIndex << " alignment: " << adata.alignments[highScoreIndex];
#endif

        // process the alignment into information that's easily usable in the vcf output (ie. breakends in reference coordinates)
#if 0
        SVCandidateAssemblyData::JumpAlignmentResultType& bestAlign(adata.alignments[adata.bestAlignmentIndex]);

        // first get each alignment associated with the correct breakend:
        Alignment* bp1AlignPtr(&bestAlign.align1);
        Alignment* bp2AlignPtr(&bestAlign.align2);


        // un-reverse alignments:
        if(isBp1Reversed)
        {
            std::reverse(bpAlignPtr)
        }

        if(isBp2AlignedFirst) std::swap(bp1AlignPtr, bp2AlignPtr);

        const pos_t align1StartPos(static_cast<pos_t>(bp1AlignPtr->alignStart));
        const pos_t align1EndPos(static_cast<pos_t>(alignEnd(&bp1AlignPtr)));

        const bool isBp1AtAlignEnd(sv.bp1.state == SVBreakendState::RIGHT_OPEN);

        adata.sv = sv;

        const pos_t bp1BreakendPos = bp1ref.get_offset() + bp1AlignPtr->;

                + adata.alignments[adata.bestAlignmentIndex].align1.
                adata.sv.bp1.interval.range.set_begin_pos()
#endif
    }
}
