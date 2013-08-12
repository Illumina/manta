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

static
bool
hasAlignedPrefix(const Alignment& al) {
	static const boost::regex re("^(\\d+)M");
	std::string cigar;
	apath_to_cigar(al.apath,cigar);
	return (boost::regex_match(cigar, re));
}

static
bool
hasAlignedSuffix(const Alignment& al) {
	static const boost::regex re("(\\d+)M$");
	std::string cigar;
	apath_to_cigar(al.apath,cigar);
	return (boost::regex_match(cigar, re));
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
