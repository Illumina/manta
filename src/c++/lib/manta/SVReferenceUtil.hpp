//
// Manta - Structural Variant and Indel Caller
// Copyright (c) 2013-2025 Illumina, Inc.
//
// This program is licensed under the terms of the Polyform strict license
//
// ***As far as the law allows, the software comes as is, without
// any warranty or condition, and the licensor will not be liable
// to you for any damages arising out of these terms or the use
// or nature of the software, under any kind of legal claim.***
//
// You should have received a copy of the PolyForm Strict License 1.0.0
// along with this program.  If not, see <https://polyformproject.org/licenses/strict/1.0.0>.
//
//

/// \file
/// \author Chris Saunders
///

#pragma once

#include "blt_util/reference_contig_segment.hpp"
#include "htsapi/bam_header_info.hpp"
#include "manta/SVCandidate.hpp"

#include <string>

/// test whether the two svCandidate breakend regions will overlap
/// after the extra reference padding has been added
bool isRefRegionOverlap(const bam_header_info& header, const pos_t extraRefEdgeSize, const SVCandidate& sv);

/// test whether the interval intersects a reference contig at all
///
bool isRefRegionValid(const bam_header_info& header, const GenomeInterval& interval);

/// Given a genome interval, attempt to add an extra buffer
/// to the interval and return the reference sequence corresponding
/// to this interval
///
/// \param[in] extraRefEdgeSize add this value to the ends of each interval prior to chromosome length
/// clipping and reference retrieval
///
/// \param[out] leadingTrim indicates how much was cut from the front of the requested interval (with edge
/// buffer)
///
/// \param[out] trailingTrim indicates how much was cut from the end of the requested interval (with edge
/// buffer)
///
void getIntervalReferenceSegment(
    const std::string&        referenceFilename,
    const bam_header_info&    header,
    const pos_t               extraRefEdgeSize,
    const GenomeInterval&     interval,
    reference_contig_segment& intervalRef,
    unsigned&                 leadingTrim,
    unsigned&                 trailingTrim);

/// alternate interface to getIntervalReferenceSegment for applications
/// where the returned trim value is not needed:
inline void getIntervalReferenceSegment(
    const std::string&        referenceFilename,
    const bam_header_info&    header,
    const pos_t               extraRefEdgeSize,
    const GenomeInterval&     interval,
    reference_contig_segment& intervalRef)
{
  unsigned leadingTrim;
  unsigned trailingTrim;
  getIntervalReferenceSegment(
      referenceFilename, header, extraRefEdgeSize, interval, intervalRef, leadingTrim, trailingTrim);
}

/// extract the reference sequence around each breakend into a reference_contig_segment
/// object
///
/// for each region, we extract the hypothetical breakend region + extraRefEdgeSize bases
/// on each side, but limit the region to [0,chrom_size-1]
///
void getSVReferenceSegments(
    const std::string&        referenceFilename,
    const bam_header_info&    header,
    const pos_t               extraRefEdgeSize,
    const SVCandidate&        sv,
    reference_contig_segment& bp1ref,
    reference_contig_segment& bp2ref,
    unsigned&                 bp1LeadingTrim,
    unsigned&                 bp1TrailingTrim,
    unsigned&                 bp2LeadingTrim,
    unsigned&                 bp2TrailingTrim);
