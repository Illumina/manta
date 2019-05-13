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
///

#pragma once

#include <cstdlib>
#include <iosfwd>
#include <string>

#include "blt_util/align_path.hpp"
#include "manta/SVBreakend.hpp"

struct SVCandidate {
  /// If false, the breakend interval is at base-pair resolution
  ///
  /// false does not mean that the interval size is zero, a precise breakend interval range represents
  /// microhomology at the breakend site
  bool isImprecise() const { return _isImprecise; }

  /// \brief Test if this SVCandidate intersects with \p rhs
  ///
  /// Two SVCandidates intersect if their breakend regions overlap in the same direction.
  ///
  /// In the schematic below, the intersecting candidate pairs are (1,2) and (2,3)
  ///
  /// Candidate1: >bp2>----------------------------------->>bp1>>
  /// Candidate2:    >>>>>>>bp1>>>>>>>--------------->>bp2>>
  /// Candidate3:               >>>bp2>>>--------------->>>>>>>bp1>>>>>>
  /// Candidate3:               <<<bp2<<<---------------<<<<<<<bp1<<<<<<
  bool isIntersect(const SVCandidate& rhs) const
  {
    return (
        (bp1.isIntersect(rhs.bp1) && bp2.isIntersect(rhs.bp2)) ||
        (bp1.isIntersect(rhs.bp2) && bp2.isIntersect(rhs.bp1)));
  }

  /// \brief Test if two SVCandidates intersect such that the two candidate's bp1 and bp2 labels match up.
  bool isIntersect1to1(const SVCandidate& rhs) const
  {
    return (bp1.isIntersect(rhs.bp1) && bp2.isIntersect(rhs.bp2));
  }

  /// \param[in] isExpandRegion If true, allow the breakpoint regions of this SVCandidate to expand to the
  /// union of this and \p rhs.
  ///
  /// \return False if the SVCandidates can't be merged because they do not intersect.
  bool merge(const SVCandidate& rhs, const bool isExpandRegion = true)
  {
    if (!isIntersect(rhs)) return false;

    if (bp1.isIntersect(rhs.bp1)) {
      bp1.merge(rhs.bp1, isExpandRegion);
      bp2.merge(rhs.bp2, isExpandRegion);
      forwardTranscriptStrandReadCount += rhs.forwardTranscriptStrandReadCount;
      reverseTranscriptStrandReadCount += rhs.reverseTranscriptStrandReadCount;
    } else {
      bp1.merge(rhs.bp2, isExpandRegion);
      bp2.merge(rhs.bp1, isExpandRegion);
      forwardTranscriptStrandReadCount += rhs.reverseTranscriptStrandReadCount;
      reverseTranscriptStrandReadCount += rhs.forwardTranscriptStrandReadCount;
    }

    _isImprecise = (isImprecise() || rhs.isImprecise());

    return true;
  }

#if 0
    void
    clear()
    {
        _isImprecise = true;
        bp1.clear();
        bp2.clear();
        insertSeq.clear();
        candidateIndex=0;
        assemblyAlignIndex=0;
        assemblySegmentIndex=0;
        isUnknownSizeInsertion = false;
        unknownSizeInsertionLeftSeq.clear();
        unknownSizeInsertionRightSeq.clear();
        forwardTranscriptStrandReadCount = 0;
        reverseTranscriptStrandReadCount = 0;
        isSingleJunctionFilter = false;
    }
#endif

  void setPrecise() { _isImprecise = false; }

  bool isForward() const { return (forwardTranscriptStrandReadCount > reverseTranscriptStrandReadCount); }

  bool isTranscriptStrandKnown() const
  {
    return (
        (std::max(forwardTranscriptStrandReadCount, reverseTranscriptStrandReadCount) + 1) /
            (std::min(forwardTranscriptStrandReadCount, reverseTranscriptStrandReadCount) + 1) >=
        2);
  }

  /// if 1 is added to the position of one breakend (within the homologous breakend range), then is 1 also
  /// added to the other breakend?
  ///
  /// if false then breakends move in opposite directions;
  bool isBreakendRangeSameShift() const { return (bp1.state != bp2.state); }

  int centerSize() const
  {
    return std::abs(bp2.interval.range.center_pos() - bp1.interval.range.center_pos());
  }

  /// Report the full spanning count
  unsigned getPostAssemblySpanningCount(const bool isRNA) const
  {
    if (!isRNA && isImprecise()) {
      return bp1.getPairCount();
    } else {
      return bp1.getSpanningCount();
    }
  }

private:
  bool _isImprecise = true;

public:
  SVBreakend bp1;
  SVBreakend bp2;

  // this is either a micro-insertion in a large-scale SV, or the inserted sequence of an actual insertion
  // in case bp1 and bp2 are on opposite strands (ie. an inversion) the insertSeq is oriented to the fwd
  // strand for bp1
  std::string insertSeq;

  /// Assembled contig sequence for outputting in VCF
  std::string contigSeq;

  // for some small indels, the alignment becomes sufficiently complex that a CIGAR string provides better
  // detail
  // (this is provided for any small SV which is more complicated than a simple insert or deletion)
  ALIGNPATH::path_t insertAlignment;

  /// Low-res candidate index number, used to generate unique SV id
  unsigned candidateIndex = 0;

  /// High-res assembly index number of alignment, used to generate unique SV id
  unsigned assemblyAlignIndex = 0;

  /// High-res assembly index number of alignment segment, used to generate unique SV id
  unsigned assemblySegmentIndex = 0;

  /// If true, the insertion hasn't been assembled all the way through
  bool isUnknownSizeInsertion = false;

  /// For an incomplete insertion, this is the known left side of the insert sequence
  std::string unknownSizeInsertionLeftSeq;

  /// For an incomplete insertion, this is the known right side of the insert sequence
  std::string unknownSizeInsertionRightSeq;

  /// Number of reads (pairs) supporting a direction from bp1 to bp2 (used for stranded RNA data)
  unsigned forwardTranscriptStrandReadCount = 0;

  /// Number of reads (pairs) directed from bp2 to bp1
  unsigned reverseTranscriptStrandReadCount = 0;

  /// If true, this sv candidate should be filtered out based on information in the candidate's own single
  /// junction.
  ///
  /// This filter may be disregarded if additional support is found for the SV candidate as part of a
  /// multi-junction event.
  bool isSingleJunctionFilter = false;
};

std::ostream& operator<<(std::ostream& os, const SVCandidate& svc);

/// \brief Specify the nature of SV evidence provided by a single DNA/RNA fragment.
///
/// A given fragment to which paired-end sequencing has been run, can provide evidence of an SV via
/// the alignment of read1 or read2 (for instance if either read directly crosses an SV breakpoint)
/// or the relative alignment of both reads, which can be used to infer that the fragment spans an
/// SV breakpoint.
///
namespace SourceOfSVEvidenceInDNAFragment {
enum index_t { UNKNOWN, READ1, READ2, READ_PAIR };

inline const char* label(const index_t i)
{
  switch (i) {
  case UNKNOWN:
    return "unknown";
  case READ1:
    return "read1";
  case READ2:
    return "read2";
  case READ_PAIR:
    return "pair";
  default:
    return "";
  }
}
}  // namespace SourceOfSVEvidenceInDNAFragment

/// \brief A specialized SVCandidate which represents an SV hypothesis generated from a single piece of
/// evidence, ie. a single SV 'observation'.
///
/// It is helpful to represent this case as a distinct specialization of SVCandidate because this allows
/// additional detail on the nature of the SV evidence to be added onto the object, (eg. "This SV candidate is
/// inferred from the poorly aligned end of the second read in a read pair").
struct SVObservation : public SVCandidate {
  SVObservation()
    : SVCandidate(),
      svEvidenceType(SVEvidenceType::UNKNOWN),
      dnaFragmentSVEvidenceSource(SourceOfSVEvidenceInDNAFragment::UNKNOWN)
  {
  }

  /// \return True if the evidence for this SV observation relies on only a single read (eg. CIGAR read
  /// alignment) or relies on both reads of a paired end observation (eg. anomalous read pair)
  bool isSingleReadSource() const
  {
    using namespace SourceOfSVEvidenceInDNAFragment;
    return ((dnaFragmentSVEvidenceSource == READ1) || (dnaFragmentSVEvidenceSource == READ2));
  }

  /// \return True if this observation is inferred from the alignment or sequence of read1
  bool isRead1Source() const
  {
    using namespace SourceOfSVEvidenceInDNAFragment;
    return (dnaFragmentSVEvidenceSource == READ1);
  }

  /// \brief The type of SV evidence (eg. anomalous read pair, large indel in CIGAR string of one read,
  /// etc....)
  SVEvidenceType::index_t svEvidenceType;
  /// \brief Ths source of the SV evidence within a single DNA fragment (eg. read1, read2, of both)
  SourceOfSVEvidenceInDNAFragment::index_t dnaFragmentSVEvidenceSource;
};

std::ostream& operator<<(std::ostream& os, const SVObservation& svc);
