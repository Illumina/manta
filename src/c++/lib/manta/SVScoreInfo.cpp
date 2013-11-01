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
/// \author Chris Saunders and Xiaoyu Chen
///

#include "manta/SVScoreInfo.hh"
#include "blt_util/log.hh"
#include "blt_util/seq_util.hh"
#include "blt_util/seq_printer.hh"

#include "boost/foreach.hpp"

#include <iostream>



SVAlignmentInfo::
SVAlignmentInfo(
    const SVCandidate& sv,
    const SVCandidateAssemblyData& assemblyData) :
    _isSpanning(assemblyData.isSpanning),
    _bp1ContigReversed(assemblyData.bporient.isBp1Reversed),
    _bp2ContigReversed(assemblyData.bporient.isBp2Reversed)
{
    // for imprecise SVs, split-read evidence won't be assigned
    if (sv.isImprecise()) return;

    contigSeq = assemblyData.extendedContigs[sv.assemblyAlignIndex];

    if (_isSpanning)
    {
        const JumpAlignmentResult<int>& alignment = assemblyData.spanningAlignments[sv.assemblyAlignIndex];

        // get offsets of breakpoints in the extended contig
        const unsigned align1Size(apath_read_length(alignment.align1.apath));
        const unsigned insertSize(alignment.jumpInsertSize);
        // the beginPos of align1 is the length of reference padding in the extended contig
        // |ref padding| + |align1| + |insert| + |align2|
        // both bp1 and bp2 include the insert and micro-homology,
        // which can avoid false split-read evidence from normal sample when the micorhomology is long

        // csaunders: leave homology size out until we're prepared downstream to
        // deal with each breakpoint as a range instead of a point value

        //unsigned homologySize = sv.bp1.interval.range.size() - 1;
        bp1ContigOffset = alignment.align1.beginPos + align1Size;
        bp2ContigOffset = alignment.align1.beginPos + align1Size + insertSize;
        if (assemblyData.bporient.isBp2AlignedFirst)
        {
            std::swap(bp1ContigOffset, bp2ContigOffset);
        }

        if (_bp1ContigReversed || _bp2ContigReversed)
        {
            revContigSeq = reverseCompCopyStr(contigSeq);
            // reset offset w.r.t. the reversed contig
            if (_bp1ContigReversed)
                bp1ContigOffset = contigSeq.size() - bp1ContigOffset - 1;
            else
                bp2ContigOffset = contigSeq.size() - bp2ContigOffset - 1;
        }

        // get reference regions
        const reference_contig_segment& bp1Ref = assemblyData.bp1ref;
        const reference_contig_segment& bp2Ref = assemblyData.bp2ref;
        bp1RefSeq = bp1Ref.seq();
        bp2RefSeq = bp2Ref.seq();
        // get offsets of breakpoints in the reference regions
        // again, both bp1 and bp2 include the micro-homology

        // csaunders: see above regarding micro-homology handling

        bp1RefOffset = sv.bp1.interval.range.begin_pos() - bp1Ref.get_offset();
        const pos_t bp2BeginPos = (sv.isBreakendRangeSameShift() ?
                                   sv.bp2.interval.range.begin_pos() :
                                   sv.bp2.interval.range.end_pos()-1);
        bp2RefOffset = bp2BeginPos - bp2Ref.get_offset();
    }
    else
    {
        // get offsets of breakpoints in the extended contig
        const AlignmentResult<int>& alignment = assemblyData.smallSVAlignments[sv.assemblyAlignIndex];
        const std::pair<unsigned, unsigned>& alignSegment = assemblyData.smallSVSegments[sv.assemblyAlignIndex][sv.assemblySegmentIndex];

        const ALIGNPATH::path_t apathTillSvStart(&alignment.align.apath[0], &alignment.align.apath[alignSegment.first]);
        const ALIGNPATH::path_t apathTillSvEnd(&alignment.align.apath[0], &alignment.align.apath[alignSegment.second+1]);

        // the beginPos of align is the length of reference padding in the extended contig
        // |ref padding| + |alignment segments|
        // both bp1 and bp2 include the insert and micro-homology,
        // which can avoid false split-read evidence from normal sample when the micorhomology is long

        // csaunders: removing microhomology

//        unsigned homologySize = sv.bp1.interval.range.size() - 1;
        bp1ContigOffset = alignment.align.beginPos + apath_read_length(apathTillSvStart);
        bp2ContigOffset = alignment.align.beginPos + apath_read_length(apathTillSvEnd);

        // get reference regions
        // only bp1ref is used for small events
        const reference_contig_segment& bp1Ref = assemblyData.bp1ref;
        bp1RefSeq = bp1Ref.seq();
        // get offsets of breakpoints in the reference regions
        // again, both bp1 and bp2 include the micro-homology

        bp1RefOffset = sv.bp1.interval.range.begin_pos() - bp1Ref.get_offset();
        bp2RefOffset = sv.bp2.interval.range.begin_pos() - bp1Ref.get_offset();
    }
}



std::ostream&
operator<<(
    std::ostream& os,
    const SVAlignmentInfo& ai)
{
    os << "Contig seq\n";
    printSeq(ai.contigSeq,os);
    os << '\n';
    os << "bp1 contig offset = " << ai.bp1ContigOffset << " bp1 contig reversed = " << ai._bp1ContigReversed << "\n";
    os << "bp2 contig offset = " << ai.bp2ContigOffset << " bp2 contig reversed = " << ai._bp2ContigReversed << "\n";
    os << "bp1RefSeq\n";
    printSeq(ai.bp1RefSeq,os);
    os << '\n';
    os << "bp2RefSeq (null for small SVs)\n";
    printSeq(ai.bp2RefSeq,os);
    os << '\n';
    os << "bp1 reference offset = " << ai.bp1RefOffset << "\n";
    os << "bp2 reference offset = " << ai.bp2RefOffset << "\n";
    return os;
}



std::ostream&
operator<<(
    std::ostream& os,
    const SVSampleAlleleInfo& sai)
{
    static const char indent('\t');
    os << "SVSampleAlleleInfo:\n"
       << indent << "bp1SpanReadCount: " << sai.bp1SpanReadCount << "\n"
       << indent << "bp2SpanReadCount: " << sai.bp2SpanReadCount << "\n"
       << indent << "spanPairCount: " << sai.spanPairCount << "\n"
       << indent << "confidentSpanningPairCount: " << sai.confidentSpanningPairCount << "\n"
       << indent << "splitReadCount: " << sai.splitReadCount << "\n"
       << indent << "confidentSplitReadCount: " << sai.confidentSplitReadCount << "\n"
       ;
    return os;
}



std::ostream&
operator<<(
    std::ostream& os,
    const SVSampleInfo& si)
{
    os << "SVSampleInfo:\n"
       << "Alt Allele " << si.alt
       << "Ref Allele " << si.ref
       ;
    return os;
}



std::ostream&
operator<<(
    std::ostream& os,
    const SVScoreInfo& ssi)
{
    os << "SVScoreInfo bp1MaxDepth=" << ssi.bp1MaxDepth << " bp2MaxDepth=" << ssi.bp2MaxDepth << "\n";
    os << "Normal sample info " << ssi.normal;
    return os;
}
