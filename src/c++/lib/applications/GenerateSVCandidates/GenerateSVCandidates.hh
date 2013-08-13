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

#pragma once

#include <string>
#include <boost/regex.hpp>

#include "manta/Program.hh"
#include "manta/SVCandidate.hh"

#include "svgraph/GenomeInterval.hh"

#include "blt_util/bam_header_util.hh"
#include "blt_util/reference_contig_segment.hh"
#include "blt_util/samtools_fasta_util.hh"
#include "blt_util/bam_header_util.hh"

#include "alignment/Alignment.hh"
#include "alignment/GlobalJumpAligner.hh"

static
void
getIntervalReferenceSegment(
    const std::string& referenceFilename,
    const bam_header_info& header,
    const pos_t extraRefEdgeSize,
    const GenomeInterval& interval,
    reference_contig_segment& intervalRef)
{
    const bam_header_info::chrom_info& chromInfo(header.chrom_data[interval.tid]);
    const std::string& chrom(chromInfo.label);
    const pos_t beginPos(std::max(0, (interval.range.begin_pos()-extraRefEdgeSize)));
    const pos_t endPos(std::min(static_cast<pos_t>(chromInfo.length), (interval.range.end_pos()+extraRefEdgeSize)));

    // get REF
    intervalRef.set_offset(beginPos);
    get_standardized_region_seq(referenceFilename,chrom,beginPos,endPos,intervalRef.seq());
}


/// extract the reference sequence around each breakend into a reference_contig_segment
/// object
///
/// for each region, we extract the hypothetical breakend region + extraRefEdgeSize bases
/// on each side
///
static
void
getSVReferenceSegments(
    const std::string& referenceFilename,
    const bam_header_info& header,
    const pos_t extraRefEdgeSize,
    const SVCandidate& sv,
    reference_contig_segment& bp1ref,
    reference_contig_segment& bp2ref)
{
    getIntervalReferenceSegment(referenceFilename,header,extraRefEdgeSize,sv.bp1.interval,bp1ref);
    getIntervalReferenceSegment(referenceFilename,header,extraRefEdgeSize,sv.bp2.interval,bp2ref);
}

// tests if prefix of aligned sequence matches target, returns length of alignment (zero if no match)
static
unsigned
hasAlignedPrefix(const Alignment& al, const unsigned minAlignContext = 0)
{
    if (al.apath.empty()) return false;
    unsigned alignLen(0);
    if (al.apath[0].type == ALIGNPATH::MATCH && al.apath[0].length >= minAlignContext)
    {
        alignLen = al.apath[0].length;
    }
    return alignLen;
}

// tests if suffix of aligned sequence matches target, returns length of alignment (zero if no match)
static
unsigned
hasAlignedSuffix(const Alignment& al, const unsigned minAlignContext = 0)
{
    if (al.apath.empty()) return false;
    size_t apLen = al.apath.size();
    unsigned alignLen(0);
    if (al.apath[apLen-1].type == ALIGNPATH::MATCH && al.apath[apLen-1].length >= minAlignContext)
    {
        alignLen = al.apath[apLen-1].length;
    }
    return alignLen;
}

static
bool
bothEndsAligned(const Alignment& al, const unsigned minAlignContext = 0)
{
    return (hasAlignedPrefix(al,minAlignContext) && hasAlignedSuffix(al,minAlignContext));
}

// check a jump alignment for consistency (only one end aligning)
static
bool
isConsistentAlignment(const JumpAlignmentResult<int>& res, const unsigned /*minAlignContext = 0*/)
{
    if (res.align1.isAligned() && res.align2.isAligned()) return true;
    return false;

    // not consistent if both unaligned
    //if (! (res.align1.isAligned() && res.align2.isAligned()) ) return false;

    // not consistent if both ends aligned for each alignment
    //if ( bothEndsAligned(res.align1) && bothEndsAligned(res.align2) ) return false;
   
	/*return ( (hasAlignedPrefix(res.align1,minAlignContext) && hasAlignedSuffix(res.align2,minAlignContext)) ||
			 (hasAlignedSuffix(res.align1,minAlignContext) && hasAlignedPrefix(res.align2,minAlignContext))
		   );*/
}

static
int
estimateBreakPointPos(const Alignment& al, const unsigned refOffset)
{
    // -1 means no breakpoint found
    int breakPointPosEstimate(-1);

    unsigned prefAlLen = hasAlignedPrefix(al);
    unsigned suffAlLen = hasAlignedSuffix(al);

    if (! (prefAlLen || suffAlLen) )
    {
        return breakPointPosEstimate;
    }

	if (prefAlLen) {
		breakPointPosEstimate = refOffset + al.alignStart + prefAlLen;
	} else if (suffAlLen) {
		breakPointPosEstimate = refOffset + alignEnd(al) - suffAlLen;
	}

    assert(breakPointPosEstimate>0);

    return breakPointPosEstimate;
}


/// generates candidate calls from graph edges
///
struct GenerateSVCandidates : public manta::Program
{

    const char*
    name() const
    {
        return "GenerateSVCandidates";
    }

    void
    runInternal(int argc, char* argv[]) const;
};
