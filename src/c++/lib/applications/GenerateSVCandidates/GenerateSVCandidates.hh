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

#include "manta/Program.hh"
#include "manta/SVCandidate.hh"

#include <string>

#include "svgraph/GenomeInterval.hh"

#include "blt_util/bam_header_util.hh"
#include "blt_util/reference_contig_segment.hh"
#include "blt_util/samtools_fasta_util.hh"
#include "blt_util/bam_header_util.hh"

#include <boost/regex.hpp>
#include "alignment/Alignment.hh"

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

// tests if prefix of aligned sequence matches target, returns length of alignment (zero if no alignment)
static
unsigned
hasAlignedPrefix(const Alignment& al, const unsigned minMatchLen) 
{
    if (al.apath.empty()) return false;
    //std::cout << "hasAlignedSuffix: " << segment_type_to_cigar_code(al.apath[0].type) << " " << al.apath[0].length << std::endl;
    unsigned alignLen(0);
    if (al.apath[0].type == ALIGNPATH::MATCH && al.apath[0].length >= minMatchLen) {
        alignLen = al.apath[0].length;
    }
    return alignLen;
}

// tests if suffix of aligned sequence matches target, returns length of alignment (zero if no alignment)
static
unsigned
hasAlignedSuffix(const Alignment& al, const unsigned minMatchLen) 
{
    if (al.apath.empty()) return false;
    size_t apLen = al.apath.size();
    //std::cout << "hasAlignedSuffix: " << segment_type_to_cigar_code(al.apath[apLen-1].type) << " " << al.apath[apLen-1].length << std::endl;
    unsigned alignLen(0);
    if (al.apath[apLen-1].type == ALIGNPATH::MATCH && al.apath[apLen-1].length >= minMatchLen) {
        alignLen = al.apath[apLen-1].length;
    }
    return alignLen;
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
