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
        const pos_t align1Size(apath_read_length(alignment.align1.apath));

        // the beginPos of align1 is the length of reference padding in the extended contig
        // |ref padding| + |align1| + |insert| + |align2|
        // both bp1 and bp2 include the insert and micro-homology,
        // which can avoid false split-read evidence from normal sample when the micorhomology is long

        // csaunders: leave homology size out until we're prepared downstream to
        // deal with each breakpoint as a range instead of a point value

        //unsigned homologySize = sv.bp1.interval.range.size() - 1;
        bp1ContigOffset = alignment.align1.beginPos + align1Size -1;
        bp2ContigOffset = bp1ContigOffset + alignment.jumpInsertSize;
        if (assemblyData.bporient.isBp2AlignedFirst)
        {
            std::swap(bp1ContigOffset, bp2ContigOffset);
        }

        if (_bp1ContigReversed || _bp2ContigReversed)
        {
            revContigSeq = reverseCompCopyStr(contigSeq);
            // reset offset w.r.t. the reversed contig
            // note this is -2 and not -1 because we're jumping to the other "side" of the breakend:
            const pos_t revSize(contigSeq.size()-2);
            if (_bp1ContigReversed)
            {
                bp1ContigOffset = revSize - bp1ContigOffset;
            }
            else
            {
                bp2ContigOffset = revSize - bp2ContigOffset;
            }
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
        bp1ContigOffset = alignment.align.beginPos + apath_read_length(apathTillSvStart) - 1;
        bp2ContigOffset = alignment.align.beginPos + apath_read_length(apathTillSvEnd) - 1;

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



bool
SVAlignmentInfo::
isMinBpEdge(
    const unsigned minEdge) const
{
    const int iminEdge(minEdge);
    if ((bp1ContigOffset+1) < iminEdge) return false;
    if ((bp2ContigOffset+1) < iminEdge) return false;
    if ((bp1RefOffset+1) < iminEdge) return false;
    if ((bp2RefOffset+1) < iminEdge) return false;

    const pos_t contigBpSize(contigSeq.size()-1);
    if ((contigBpSize - bp1ContigOffset) < iminEdge) return false;
    if ((contigBpSize - bp2ContigOffset) < iminEdge) return false;

    const pos_t bp1RefSize(bp1ReferenceSeq().size());
    if ((bp1RefSize - 1 - bp1RefOffset) < iminEdge) return false;

    const pos_t bp2RefSize(bp2ReferenceSeq().size());
    if ((bp2RefSize - 1 - bp2RefOffset) < iminEdge) return false;

    return true;
}



static
void
dumpSeq(
    const char* label,
    const std::string& seq,
    std::ostream& os)
{
    os << label << " size/seq: " << seq.size() << '\n';
    printSeq(seq,os);
    os << '\n';
}



std::ostream&
operator<<(
    std::ostream& os,
    const SVAlignmentInfo& ai)
{
    os << "SVAlignmentInfo: isSpanning: " << ai.isSpanning() << '\n';
    dumpSeq("Contig",ai.contigSeq,os);
    dumpSeq("Rev Contig",ai.revContigSeq,os);
    os << "bp1 contig offset = " << ai.bp1ContigOffset << " bp1 contig reversed = " << ai._bp1ContigReversed << '\n';
    os << "bp2 contig offset = " << ai.bp2ContigOffset << " bp2 contig reversed = " << ai._bp2ContigReversed << '\n';
    dumpSeq("bp1RefSeq",ai.bp1RefSeq,os);
    if (ai.isSpanning())
    {
        dumpSeq("bp2RefSeq",ai.bp2RefSeq,os);
    }
    os << "bp1 reference offset = " << ai.bp1RefOffset << '\n';
    os << "bp2 reference offset = " << ai.bp2RefOffset << '\n';
    return os;
}



std::ostream&
operator<<(
    std::ostream& os,
    const SVSampleAlleleInfo& sai)
{
    static const char indent('\t');
    os << "SVSampleAlleleInfo:\n"
       << indent << "confidentSpanningPairCount: " << sai.confidentSpanningPairCount << '\n'
       << indent << "confidentSemiMappedSpanningPairCount: " << sai.confidentSemiMappedSpanningPairCount << '\n'
       << indent << "splitReadCount: " << sai.splitReadCount << '\n'
       << indent << "confidentSplitReadCount: " << sai.confidentSplitReadCount << '\n'
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
    os << "SVScoreInfo bp1MaxDepth=" << ssi.bp1MaxDepth << " bp2MaxDepth=" << ssi.bp2MaxDepth << '\n'
       << "SVScoreInfo bp1MQ0Frac=" << ssi.bp1MQ0Frac << " bp2MQ0Frac=" << ssi.bp2MQ0Frac << '\n'
       << "Normal sample info " << ssi.normal;
    return os;
}

