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



#if 0
static
void
trimOverlappingRange(
    known_pos_range2& rA,
    known_pos_range2& rB)
{
    // put ranges in order:
    known_pos_range2* r1(&rA);
    known_pos_range2* r2(&rB);

    if (r1->begin_pos() > r2->begin_pos()) std::swap(r1, r2);

    const pos_t overlap(r1->end_pos()-r2->begin_pos());
    if (overlap <= 0) return;

    r1->set_end_pos(r1->end_pos()-(overlap/2));
    r2->set_begin_pos(r1->end_pos());
}
#endif



/// produce the reference extraction interval only
static
void
getBpReferenceInterval(
    const bam_header_info& header,
    const pos_t extraRefEdgeSize,
    const GenomeInterval& bpInterval,
    GenomeInterval& refInterval)
{
    const bam_header_info::chrom_info& chromInfo(header.chrom_data[bpInterval.tid]);

    const pos_t beginPos(std::max(0, (bpInterval.range.begin_pos()-extraRefEdgeSize)));
    const pos_t endPos(std::min(static_cast<pos_t>(chromInfo.length), (bpInterval.range.end_pos()+extraRefEdgeSize)));

    refInterval.tid = bpInterval.tid;
    refInterval.range.set_begin_pos(beginPos);
    refInterval.range.set_end_pos(endPos);
}


/// given a reference extraction interval, produce the corresponding ref contig segment
static
void
getIntervalReferenceSegment(
    const std::string& referenceFilename,
    const bam_header_info& header,
    const GenomeInterval& refInterval,
    reference_contig_segment& intervalRefSeq)
{
    const bam_header_info::chrom_info& chromInfo(header.chrom_data[refInterval.tid]);
    const std::string& chrom(chromInfo.label);

    // get REF
    const known_pos_range2& range(refInterval.range);
    intervalRefSeq.set_offset(range.begin_pos());

    // note: begin and end pos follow Manta's closed-open bpInterval conventions (a la bedtools,
    // but the ref function below takes closed-closed endpoints, so we subract one from endPos
    get_standardized_region_seq(referenceFilename, chrom, range.begin_pos(), (range.end_pos()-1), intervalRefSeq.seq());

    assert(static_cast<pos_t>(intervalRefSeq.seq().size()) == (static_cast<pos_t>(range.size())));
}



void
getIntervalReferenceSegment(
    const std::string& referenceFilename,
    const bam_header_info& header,
    const pos_t extraRefEdgeSize,
    const GenomeInterval& bpInterval,
    reference_contig_segment& intervalRefSeq)
{
    GenomeInterval refInterval;
    getBpReferenceInterval(header, extraRefEdgeSize, bpInterval, refInterval);

    getIntervalReferenceSegment(referenceFilename, header, refInterval, intervalRefSeq);
}



bool
isRefRegionOverlap(
    const bam_header_info& header,
    const pos_t extraRefEdgeSize,
    const SVCandidate& sv)
{
    if (sv.bp1.interval.tid != sv.bp2.interval.tid) return false;
    GenomeInterval bp1RefInterval;
    GenomeInterval bp2RefInterval;
    getBpReferenceInterval(header,extraRefEdgeSize,sv.bp1.interval,bp1RefInterval);
    getBpReferenceInterval(header,extraRefEdgeSize,sv.bp2.interval,bp2RefInterval);

    return (bp1RefInterval.isIntersect(bp2RefInterval));
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
    GenomeInterval bp1RefInterval;
    GenomeInterval bp2RefInterval;
    getBpReferenceInterval(header,extraRefEdgeSize,sv.bp1.interval,bp1RefInterval);
    getBpReferenceInterval(header,extraRefEdgeSize,sv.bp2.interval,bp2RefInterval);

    // allow overlap (best performance in case of breakends in opposite orientations...:
#if 0
    // check that the two reference regions do not overlap
    if (bp1RefInterval.isIntersect(bp2RefInterval))
    {
        // rare case, trim intervals so that they become non-overlapping:
        trimOverlappingRange(bp1RefInterval.range, bp2RefInterval.range);
    }
#endif

    getIntervalReferenceSegment(referenceFilename, header, bp1RefInterval, bp1ref);
    getIntervalReferenceSegment(referenceFilename, header, bp2RefInterval, bp2ref);
}
