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

#include "boost/foreach.hpp"

#include <iostream>



SVAlignmentInfo::
SVAlignmentInfo(
    const SVCandidate& sv,
    const SVCandidateAssemblyData& assemblyData)
{
    // consider 2-locus events first
    // TODO: to add local assembly later

    // for imprecise SVs, split-read evidence won't be assigned
    if ((assemblyData.isSpanning) &&
        (!sv.isImprecise()))
    {
        contigSeq = assemblyData.extendedContigs[assemblyData.bestAlignmentIndex];
        const JumpAlignmentResult<int>& alignment = assemblyData.spanningAlignments[assemblyData.bestAlignmentIndex];

        // get offsets of breakpoints in the contig
        const unsigned align1Size(apath_read_length(alignment.align1.apath));
        const unsigned insertSize(alignment.jumpInsertSize);
        // the beginPos of align1 is the length of reference padding in the extended contig
        // |ref padding| + |align1| + |insert| + |align2|
        bp1ContigOffset = alignment.align1.beginPos + align1Size - 1;
        bp2ContigOffset = alignment.align1.beginPos + align1Size + insertSize;
        if (assemblyData.isBp2AlignedFirst)
        {
            std::swap(bp1ContigOffset, bp2ContigOffset);
        }

        bp1ContigReversed = assemblyData.isBp1Reversed;
        bp2ContigReversed = assemblyData.isBp2Reversed;

        if (bp1ContigReversed || bp2ContigReversed)
        {
            revContigSeq = reverseCompCopyStr(contigSeq);
            // reset offset w.r.t. the reversed contig
            if (bp1ContigReversed)
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
        bp1RefOffset = sv.bp1.interval.range.begin_pos() - bp1Ref.get_offset();
        const pos_t bp2BeginPos = (sv.isBreakendRangeSameShift() ?
                                   sv.bp2.interval.range.begin_pos() :
                                   sv.bp2.interval.range.end_pos()-1);
        bp2RefOffset = bp2BeginPos - bp2Ref.get_offset();
    }
}



std::ostream&
operator<<(
    std::ostream& os,
    const SVAlignmentInfo& ai)
{
    os << "Contig seq\n" << ai.contigSeq << "\n";
    os << "bp1 contig offset = " << ai.bp1ContigOffset << " bp1 contig reversed = " << ai.bp1ContigReversed << "\n";
    os << "bp2 contig offset = " << ai.bp2ContigOffset << " bp2 contig reversed = " << ai.bp2ContigReversed << "\n";
    os << "bp1RefSeq\n" << ai.bp1RefSeq << "\n";
    os << "bp2RefSeq (null for small SVs)\n" << ai.bp2RefSeq << "\n";
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
       << indent << "splitReadCount: " << sai.splitReadCount << "\n"
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
    BOOST_FOREACH(const std::string& filter, ssi.filters)
    {
        os << " " << filter;
    }
    return os;
}
