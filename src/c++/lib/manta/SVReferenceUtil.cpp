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

#include "manta/SVReferenceUtil.hh"

#include "blt_util/samtools_fasta_util.hh"



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
