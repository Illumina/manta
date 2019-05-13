//
// Manta - Structural Variant and Indel Caller
// Copyright (c) 2013-2019 Illumina, Inc.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//

/// \file
/// \author Chris Saunders
/// \author Ole Schulz-Trieglaff
/// \author Bret Barnes
/// \author Felix Schlesinger
///

#pragma once

#include "blt_util/SimpleAlignment.hpp"
#include "blt_util/blt_types.hpp"
#include "blt_util/reference_contig_segment.hpp"
#include "htsapi/bam_record.hpp"

/// \brief Searches for reads with poorly aligned ends indicative of a possible SV or indel breakpoint.
///
/// The input read is examined for evidence of poor alignment on both the leading and trailing ends.
/// If poor edge alignment is found, the read is refered to as "semi-aligned' because in this state
/// part of the read is well aligned while one or both edges of the read are poorly aligned. A semi-aligned
/// read is distinct from a read which is poorly aligned throughout.
///
/// To do so the following steps are taken:
/// 1. Reads where there is suspected adaptor 'read-through' due to short fragments size are filtered out.
/// 2. Reads which overlap their mate read are optionally filtered out (per function argument \p
/// useOverlapPairEvidence).
/// 3. Any soft-clipped segments in the input reed are 'unrolled' to the match state.
/// 4. The length of poorly aligned sequence on each read edge is found.
/// 5. Cases where the entire read is poorly aligned are filtered out.
/// 6. For overlapping read pairs, any poorly aligned edges on the "interior" of the DNA fragment are filtered
/// out. (see elaboration on this point below)
/// 7. Poorly aligned edges are filtered out if they do not have high basecall quality, at least \p
/// minHighBasecallQualityFraction of the edge basecalls must have basecall quality >= \p minBasecallQuality.
///
/// The logic behind step (6) above is as follows. Consider an overlapping read pair represented schematically
/// as:
//
/// read1 |--------------->
/// read2         <--------------|
///
/// Considering this DNA fragment as a whole, we expect a poor alignment on either the left end of read1 or
/// the right end of read2 to indicate that a breakend opens to the left or right side of the fragment,
/// respectively.
///
/// In the following example we add "S" to signify a portion of the read which is poorly aligned:
//
/// read1 |------------SSS>
/// read2         <----SSSSSSSSSS|
///
/// Now the DNA fragment is showing clear evidence of a breakend opening to the right side of the fragment. If
/// we considered the semi-aligned sequence on the right side or read1 as breakend evidence, we would fail to
/// recognize that this observation is not independent of the similar (but better) evidence from read2, so in
/// this case, step 6 would filter out the evidence from read1, but keep the evidence from read2.
///
/// \param[in] bamRead The input read being evaluated.
///
/// \param[in] bamAlign A simplified representation of the bamRead's alignment. This information can be
/// derived from bamRead but is provided as a (presumed) cache optimization.
///
/// \param[in] refSeq Local reference sequence
///
/// \param[in] useOverlapPairEvidence When false, filter out read pairs which overlap.
///
/// \param[in] isAggressiveAdaptorCheck When true, filter out read pairs that might run into adapters based on
/// an aggressive check
///
/// \param[out] leadingEdgePoorAlignmentLength Length of possible breakend-associated poor alignment on the
/// read's leading edge.
///
/// \param[out] leadingEdgeRefPos The reference position of the first aligned base after the leading poorly
/// aligned segment.
///
/// \param[out] trailingEdgePoorAlignmentLength Length of possible breakend-associated poor alignment on the
/// read's trailing edge.
///
/// \param[out] trailingEdgeRefPos The reference position of the last aligned base before the trailing poorly
/// aligned segment.
///
/// \param[in] minBasecallQuality See \p minHighBasecallQualityFraction
///
/// \param[in] minHighBasecallQualityFraction Poorly aligned edge must have \p minHighBasecallQualityFraction
/// fraction of edge basecalls with basecall quality of \p minBasecallQuality or greater, otherwise the edge
/// is filtered out.
///
void getSVBreakendCandidateSemiAligned(
    const bam_record&               bamRead,
    const SimpleAlignment&          bamAlign,
    const reference_contig_segment& refSeq,
    const bool                      useOverlapPairEvidence,
    unsigned&                       leadingEdgePoorAlignmentLength,
    pos_t&                          leadingEdgeRefPos,
    unsigned&                       trailingEdgePoorAlignmentLength,
    pos_t&                          trailingEdgeRefPos,
    const uint8_t                   minBasecallQuality             = 20,
    const float                     minHighBasecallQualityFraction = 0.75);

/// \brief Simplified interface to ::getSVBreakendCandidateSemiAligned
///
/// This version presents a subset of the parameters used by the full interface above.
inline void getSVBreakendCandidateSemiAlignedSimple(
    const bam_record&               bamRead,
    const SimpleAlignment&          bamAlign,
    const reference_contig_segment& refSeq,
    const bool                      useOverlapPairEvidence,
    unsigned&                       leadingMismatchLen,
    unsigned&                       trailingMismatchLen)
{
  pos_t leadingRefPos(0), trailingRefPos(0);
  getSVBreakendCandidateSemiAligned(
      bamRead,
      bamAlign,
      refSeq,
      useOverlapPairEvidence,
      leadingMismatchLen,
      leadingRefPos,
      trailingMismatchLen,
      trailingRefPos);
}
