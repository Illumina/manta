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

#include "SVLocusScannerSemiAligned.hpp"
#include "htsapi/bam_record_util.hpp"

#include "boost/foreach.hpp"

//#define DEBUG_SEMI_ALIGNED

#ifdef DEBUG_SEMI_ALIGNED
#include <iostream>
#include "blt_util/log.hpp"
#endif

/// \brief Base matching criteria for the purpose of finding poorly aligned read edge lengths.
///
/// For this application we treat all ambiguous bases as matching. Note that, more generally, we usually want
/// ambiguous basecalls to be counted as mismatches, but for the poorly-aligned metric, we don't want N's to
/// nominate read segments as SV-associated because of a string of Ns, so we treat these as matches.
///
static bool isBaseMatchForPoorAlignmentTest(const char a, const char b)
{
  if ((a == 'N') || (b == 'N')) return true;
  return (a == b);
}

/// \brief Leading edge helper function for ::edgePoorAlignmentLength.
static void leadingEdgePoorAlignmentLength(
    const SimpleAlignment&          bamAlign,
    const bam_seq&                  querySeq,
    const reference_contig_segment& refSeq,
    const unsigned                  contiguousMatchCount,
    unsigned&                       leadingLength,
    pos_t&                          leadingRefPos)
{
  using namespace ALIGNPATH;

  assert(contiguousMatchCount != 0);

  pos_t readIndex(0);
  pos_t refIndex(bamAlign.pos);

  leadingLength = 0;
  leadingRefPos = refIndex;

  unsigned matchLength(0);
  for (const path_segment& ps : bamAlign.path) {
    if (is_segment_align_match(ps.type)) {
      for (unsigned segPos(0); segPos < ps.length; ++segPos) {
        if (isBaseMatchForPoorAlignmentTest(
                querySeq.get_char(readIndex + segPos), refSeq.get_base(refIndex + segPos))) {
          matchLength++;

          if (matchLength >= contiguousMatchCount) {
            leadingLength = (readIndex + segPos) - (matchLength - 1);
            leadingRefPos = (refIndex + segPos) - (matchLength - 1);
            return;
          }
        } else {
          matchLength = 0;
        }
      }
    } else if (is_segment_type_indel(ps.type)) {
      matchLength = 0;
    }

    if (is_segment_type_read_length(ps.type)) readIndex += ps.length;
    if (is_segment_type_ref_length(ps.type)) refIndex += ps.length;
  }

  leadingLength = readIndex;
  leadingRefPos = refIndex;
}

/// \brief Trailing edge helper function for ::edgePoorAlignmentLength
static void trailingEdgePoorAlignmentLength(
    const SimpleAlignment&          bamAlign,
    const bam_seq&                  querySeq,
    const reference_contig_segment& refSeq,
    const unsigned                  contiguousMatchCount,
    unsigned&                       trailingLength,
    pos_t&                          trailingRefPos)
{
  using namespace ALIGNPATH;

  assert(contiguousMatchCount != 0);

  const pos_t readSize(querySeq.size());

  pos_t readIndex(readSize - 1);
  pos_t refIndex(bamAlign.pos + apath_ref_length(bamAlign.path) - 1);

  unsigned matchLength(0);
  BOOST_REVERSE_FOREACH(const path_segment& ps, bamAlign.path)
  {
    if (is_segment_align_match(ps.type)) {
      for (unsigned segPos(0); segPos < ps.length; ++segPos) {
        if (isBaseMatchForPoorAlignmentTest(
                querySeq.get_char(readIndex - segPos), refSeq.get_base(refIndex - segPos))) {
          matchLength++;

          if (matchLength >= contiguousMatchCount) {
            trailingLength = (readSize - (readIndex - segPos)) - matchLength;
            trailingRefPos = (refIndex - segPos) + matchLength;
            return;
          }
        } else {
          matchLength = 0;
        }
      }
    } else if (is_segment_type_indel(ps.type)) {
      matchLength = 0;
    }

    if (is_segment_type_read_length(ps.type)) readIndex -= ps.length;
    if (is_segment_type_ref_length(ps.type)) refIndex -= ps.length;
  }

  trailingLength = readSize - (readIndex + 1);
  trailingRefPos = refIndex + 1;
}

/// \brief Report the number of bases from each edge of the read which are poorly aligned.
///
/// The number of poorly aligned bases are determined by finding the number of bases which must be traversed
/// before encountering a contiguous matching segment of specified length. A contiguous matching segment is
/// a continuous block of \p contiguousMatchCount basecall matches between the read and the reference
/// sequence.
///
/// Example:
/// ```
/// Given a contiguousMatchCount of 3, and
/// Query: ABAAAAAABBAAB
///   Ref; AAAAAAAAAAAAA
///
/// leadingLength would be 2 and trailingLength would be 5
/// ```
///
/// \param[in] contiguousMatchCount Minimum length of contiguous matching block used to define edge mismatch
/// length.
///
/// \param[out] leadingLength Length in read coordinates from the start of the read to the base preceding the
/// first matching segment.
///
/// \param[out] leadingRefPos Reference position of read start after removing leadingLength from the beginning
/// of the read.
///
/// \param[out] trailingLength Length in read coordinates from the end of the read to the base following the
/// last matching segment.
///
/// \param[out] trailingRefPos Reference position of read end after removing trailingLength from the end of
/// the read.
///
static void edgePoorAlignmentLength(
    const SimpleAlignment&          bamAlign,
    const bam_seq&                  querySeq,
    const reference_contig_segment& refSeq,
    const unsigned                  contiguousMatchCount,
    unsigned&                       leadingLength,
    pos_t&                          leadingRefPos,
    unsigned&                       trailingLength,
    pos_t&                          trailingRefPos)
{
  leadingEdgePoorAlignmentLength(
      bamAlign, querySeq, refSeq, contiguousMatchCount, leadingLength, leadingRefPos);
  trailingEdgePoorAlignmentLength(
      bamAlign, querySeq, refSeq, contiguousMatchCount, trailingLength, trailingRefPos);

  // Verify that neither edge's poorly aligned length exceeds the read length.
  const unsigned readSize(querySeq.size());
  assert(leadingLength <= readSize);
  assert(trailingLength <= readSize);
}

void getSVBreakendCandidateSemiAligned(
    const bam_record&               bamRead,
    const SimpleAlignment&          bamAlign,
    const reference_contig_segment& refSeq,
    const bool                      useOverlapPairEvidence,
    unsigned&                       leadingEdgePoorAlignmentLength,
    pos_t&                          leadingEdgeRefPos,
    unsigned&                       trailingEdgePoorAlignmentLength,
    pos_t&                          trailingEdgeRefPos,
    const uint8_t                   minBasecallQuality,
    const float                     minHighBasecallQualityFraction)
{
  static const unsigned contiguousMatchCount(5);

  leadingEdgePoorAlignmentLength  = 0;
  leadingEdgeRefPos               = 0;
  trailingEdgePoorAlignmentLength = 0;
  trailingEdgeRefPos              = 0;

  const bool isOverlappingReadPair(is_overlapping_pair(bamRead, bamAlign));
  if (isOverlappingReadPair) {
    if ((!useOverlapPairEvidence) || (is_adapter_pair(bamRead))) return;
  }

  using namespace ALIGNPATH;
  const bam_seq querySeq(bamRead.get_bam_read());

  const uint8_t* qual(bamRead.qual());
  const unsigned readSize(bamRead.read_size());

  // Create a new alignment with all soft-clip sections unrolled to a matched alignment state.
  const SimpleAlignment matchedAlignment(matchifyEdgeSoftClip(bamAlign));

  // Get the poorly aligned length for the leading and trailing edge of the input read
  unsigned leadingEdgePoorAlignmentLengthTmp(0);
  unsigned trailingEdgePoorAlignmentLengthTmp(0);
  edgePoorAlignmentLength(
      matchedAlignment,
      querySeq,
      refSeq,
      contiguousMatchCount,
      leadingEdgePoorAlignmentLengthTmp,
      leadingEdgeRefPos,
      trailingEdgePoorAlignmentLengthTmp,
      trailingEdgeRefPos);

  // Filter out cases where the entire read is poorly aligned
  if ((leadingEdgePoorAlignmentLengthTmp + trailingEdgePoorAlignmentLengthTmp) >= readSize) return;

  // Only report the leading or trailing poorly aligned edge length for edges with
  // high basecall qualities.
  //
  // Additionally, if read pair is overlapping, don't report the edges in the interior of the
  // read fragment as poorly aligned.

  if (0 != leadingEdgePoorAlignmentLengthTmp) {
    if ((!isOverlappingReadPair) || bamRead.isSASplit() || bamRead.is_fwd_strand()) {
      unsigned highQualityCount(0);
      for (unsigned pos(0); pos < leadingEdgePoorAlignmentLengthTmp; ++pos) {
        if (qual[pos] >= minBasecallQuality) ++highQualityCount;
      }
      if ((static_cast<float>(highQualityCount) / (leadingEdgePoorAlignmentLengthTmp)) >=
          minHighBasecallQualityFraction) {
        leadingEdgePoorAlignmentLength = leadingEdgePoorAlignmentLengthTmp;
      }
    }
#ifdef DEBUG_SEMI_ALIGNED
    else {
      log_os << " Overlapping_pair leading"
             << " read qname=" << bamRead.qname() << "\n";
    }
#endif
  }

  if (0 != trailingEdgePoorAlignmentLengthTmp) {
    if ((!isOverlappingReadPair) || bamRead.isSASplit() || (!bamRead.is_fwd_strand())) {
      unsigned highQualityCount(0);
      for (unsigned pos(0); pos < trailingEdgePoorAlignmentLengthTmp; ++pos) {
        if (qual[readSize - pos - 1] >= minBasecallQuality) ++highQualityCount;
      }
      if ((static_cast<float>(highQualityCount) / (trailingEdgePoorAlignmentLengthTmp)) >=
          minHighBasecallQualityFraction) {
        trailingEdgePoorAlignmentLength = trailingEdgePoorAlignmentLengthTmp;
      }
    }
#ifdef DEBUG_SEMI_ALIGNED
    else {
      log_os << "Overlapping_pair trailing"
             << " read qname=" << bamRead.qname() << "\n";
    }
#endif
  }
}
