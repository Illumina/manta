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

    // assemble contig spanning the breakend:
    _svAssembler.assembleSVBreakends(sv.bp1, sv.bp2, adata.contigs);


    // how much additional reference sequence should we extract from around
    // each side of the breakend region?
    static const pos_t extraRefEdgeSize(300);

    // min alignment context
    //const unsigned minAlignContext(4);
    // don't align contigs shorter than this
    static const unsigned minContigLen(75);

    reference_contig_segment bp1ref,bp2ref;
    getSVReferenceSegments(_opt.referenceFilename, _header, extraRefEdgeSize, sv, bp1ref, bp2ref);
    const std::string& bp1RefStr(bp1ref.seq());
    const std::string& bp2RefStr(bp2ref.seq());

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

        // QC contigs prior to alignment:
        if (contig.seq.size() < minContigLen) continue;

        JumpAlignmentResult<int> alignment(adata.alignments[contigIndex]);

        _aligner.align(contig.seq.begin(),contig.seq.end(),
                      bp1RefStr.begin(),bp1RefStr.end(),
                      bp2RefStr.begin(),bp2RefStr.end(),
                      alignment);

#ifdef DEBUG_REFINER
        log_os << "cid: " << contigIndex << " alignment: " << alignment;
#endif

        // check which jump alignment has anchors/alignments to both breakends
        bool isAlignmentGood(alignment.align1.isAligned() && alignment.align2.isAligned());

        if(! isAlignmentGood) continue;

        if((! isHighScore) || (alignment.score > adata.alignments[highScoreIndex].score))
        {
            highScoreIndex=contigIndex;
        }
    }

    // set any additional QC steps before deciding an alignment is usable:


    // ok, passed QC -- mark the high-scoring alignment as usable:
    if(isHighScore)
    {
        adata.isBestAlignment = true;
        adata.bestAlignmnetIndex = highScoreIndex;
#ifdef DEBUG_REFINER
        log_os << "highscoreid: " << highScoreIndex << " alignment: " << adata.alignments[highScoreIndex];
#endif
    }
}
