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
/// \author Chris Saunders
///

#include "SVCandidateAssemblyRefiner.hh"

#include "alignment/AlignmentUtil.hh"
#include "blt_util/samtools_fasta_util.hh"
#include "manta/SVLocusAssembler.hh"
#include "manta/SVReferenceUtil.hh"



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
    //log_os << "sv: " << sv;

    adata.clear();

    // assemble contig spanning the breakend:
    _svAssembler.assembleSVBreakends(sv.bp1, sv.bp2, adata.contigs);
#if 1

    // how much additional reference sequence should we extract from around
    // each side of the breakend region?
    static const pos_t extraRefEdgeSize(300);

    // min alignment context
    //const unsigned minAlignContext(4);
    // don't align contigs shorter than this
    const unsigned minContigLen(75);

    reference_contig_segment bp1ref,bp2ref;
    getSVReferenceSegments(_opt.referenceFilename, _header, extraRefEdgeSize, sv, bp1ref, bp2ref);
    const std::string bp1RefStr(bp1ref.seq());
    const std::string bp2RefStr(bp2ref.seq());

    BOOST_FOREACH(const AssembledContig& contig, adata.contigs)
    {
        if (contig.seq.size() < minContigLen) continue;

        // JumpAligner allows only one jump per alignment, so we need to align
        // two times:
        JumpAlignmentResult<int> forwardRes;
        JumpAlignmentResult<int> backwardRes;

        _aligner.align(contig.seq.begin(),contig.seq.end(),
                      bp1RefStr.begin(),bp1RefStr.end(),
                      bp2RefStr.begin(),bp2RefStr.end(),
                      forwardRes);

        _aligner.align(contig.seq.begin(),contig.seq.end(),
                      bp2RefStr.begin(),bp2RefStr.end(),
                      bp1RefStr.begin(),bp1RefStr.end(),
                      backwardRes);

        // check which jump alignment has anchors/alignments to both breakends
        bool forwardGood(forwardRes.align1.isAligned() && forwardRes.align2.isAligned());
        bool backwardGood(backwardRes.align1.isAligned() && backwardRes.align2.isAligned());

        if (forwardGood && !backwardGood)
        {
            adata.alignments.push_back(forwardRes);
        }
        else if (!forwardGood && backwardGood)
        {
            adata.alignments.push_back(backwardRes);
        }
        else
        {
            // tie, store jump alignment with higher score
            adata.alignments.push_back(forwardRes.score>=backwardRes.score ? forwardRes : backwardRes);
        }
    }
#endif
}
