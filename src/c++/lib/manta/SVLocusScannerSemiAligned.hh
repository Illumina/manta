// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Manta
// Copyright (c) 2013-2015 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/sequencing/licenses/>
//

///
/// \author Chris Saunders
/// \author Ole Schulz-Trieglaff
/// \author Bret Barnes
///

#pragma once

#include "blt_util/blt_types.hh"
#include "blt_util/reference_contig_segment.hh"
#include "blt_util/SimpleAlignment.hh"
#include "htsapi/bam_record_util.hh"


/// searches for poorly aligned read ends
///
/// search is based on high-quality mismatches, any soft-clipped sections
/// will be realigned before searching for semi-aligned sections
///
/// \param[in] minQ
/// \param[in] minQFrac this fraction of bases must have qual>=minQ within the clipped region
///
void
getSVBreakendCandidateSemiAligned(
    const bam_record& bamRead,
    const SimpleAlignment& bamAlign,
    const reference_contig_segment& refSeq,
    unsigned& leadingMismatchLen,
    pos_t& leadingRefPos,
    unsigned& trailingMismatchLen,
    pos_t& trailingRefPos,
    const uint8_t minQ = 20,
    const float minQFrac = 0.75);


/// simplified interface to the full function above:
inline
void
getSVBreakendCandidateSemiAlignedSimple(
    const bam_record& bamRead,
    const SimpleAlignment& bamAlign,
    const reference_contig_segment& refSeq,
    unsigned& leadingMismatchLen,
    unsigned& trailingMismatchLen)
{
    pos_t leadingRefPos(0), trailingRefPos(0);
    getSVBreakendCandidateSemiAligned(
        bamRead, bamAlign, refSeq,
        leadingMismatchLen, leadingRefPos,
        trailingMismatchLen, trailingRefPos);
}
