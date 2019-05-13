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

#pragma once

#include <iosfwd>
#include <string>
#include <unordered_map>
#include <vector>

#include "alignment/GlobalAligner.hpp"
#include "alignment/GlobalJumpAligner.hpp"
#include "assembly/AssembledContig.hpp"
#include "blt_util/reference_contig_segment.hpp"
#include "manta/SVCandidate.hpp"

/// \brief Minimum set of information required to describe bp transformations between SVCandidate and its
/// corresponding contig alignment
///
struct BPOrientation {
  void clear()
  {
    isBp2AlignedFirst       = false;
    isBp1Reversed           = false;
    isBp2Reversed           = false;
    isBp1First              = true;
    isTranscriptStrandKnown = false;
  }

  /// should the contig on the fwd strand align bp2->bp1 (true) or bp1->bp2 (false)?
  bool isBp2AlignedFirst = false;
  bool isBp1Reversed     = false;  ///< should all bp1 reads be reversed for the contig to assemble correctly?
  bool isBp2Reversed     = false;  ///< should all bp2 reads be reversed for the contig to assemble correctly?

  /// Is this candidate oriented from bp1 to bp2 (used in RNA)? Valid if isTranscriptStrandKnown is true
  bool isBp1First = true;

  bool isTranscriptStrandKnown = false;  ///< Do we know the strand for this candidate (RNA)
};

struct LargeInsertionInfo {
  void clear()
  {
    isLeftCandidate  = false;
    isRightCandidate = false;
    contigOffset     = 0;
    refOffset        = 0;
    score            = 0;
  }

  bool     isLeftCandidate  = false;
  bool     isRightCandidate = false;
  unsigned contigOffset     = 0;  ///< If candidate, how far into the contig is the breakend?

  /// If candidate, how far from the start of the contig alignment is the breakend on reference?
  unsigned refOffset = 0;

  int score = 0;  ///< What is the alignment score of the contig up to the insertion breakpoint?
};

std::ostream& operator<<(std::ostream& os, const LargeInsertionInfo& lii);

struct RemoteReadPayload {
  RemoteReadPayload() : readNo(0) {}

  RemoteReadPayload(const int initReadNo, const std::string& initReadSeq)
    : readNo(initReadNo), readSeq(initReadSeq)
  {
  }

  uint8_t     readNo;  ///< Read no of the remote read, ie. the readno matching readSeq
  std::string readSeq;
};

typedef std::unordered_map<std::string, RemoteReadPayload> RemoteReadCache;

/// \brief Assembly data pertaining to a specific SV candidate
///
/// Assembly starts from a low-resolution SV candidate. This holds
/// any persistent data related to the assembly process, such as data
/// useful for scoring the candidate.
///
/// As a future design note -- it may be possible that the candidate is broken
/// into multiple candidates during assembly-based hypothesis refinement, so
/// this struct may cover multiple refined candidates (but always only one input
/// candidate.
///
/// Also note that this class is representing both large scale/spanning SV's using the locus 'jump' aligner
/// and small-scale local assemblies. This should probably be refactored into two parts, but
/// it's hard to see the right strategy/interface until the scoring modules reach greater maturity
/// (the scoring modules are the primary non-local consumer of information from this struct)
///
struct SVCandidateAssemblyData {
  void clear()
  {
    contigs.clear();
    isCandidateSpanning = false;
    isSpanning          = false;
    bporient.clear();
    extendedContigs.clear();
    smallSVAlignments.clear();
    spanningAlignments.clear();
    smallSVSegments.clear();
    largeInsertInfo.clear();
    remoteReads.clear();
    bestAlignmentIndex = 0;
    bp1ref.clear();
    bp2ref.clear();
    svs.clear();
    isOverlapSkip = false;
  }

  typedef AlignmentResult<int>     SmallAlignmentResultType;
  typedef JumpAlignmentResult<int> JumpAlignmentResultType;

  // Here 'segments' only pertain to small-event alignments:
  typedef std::pair<unsigned, unsigned>     CandidateSegmentType;
  typedef std::vector<CandidateSegmentType> CandidateSegmentSetType;

  Assembly contigs;  ///< assembled contigs for both breakpoints

  /// \brief If true, then before assembly this was a 2-locus event (spanning) instead of a local-assembly
  ///
  /// Note that isCandidateSpanning can be != isSpanning in cases where the breakends are so close that the
  /// assembler decides to fuse them and treat a spanning hypothesis as non-spanning:
  bool isCandidateSpanning = false;

  /// If true, then at assembly time this was treated as a 2-locus event (spanning) instead of a
  /// local-assembly
  bool isSpanning = false;

  BPOrientation bporient;

  /// Assembled contigs with additional reference sequence padding added onto each edge
  std::vector<std::string> extendedContigs;

  /// Contig smallSV alignments, one per contig, may be empty
  std::vector<SmallAlignmentResultType> smallSVAlignments;

  /// Contig spanning alignments, one per contig, may be empty
  std::vector<JumpAlignmentResultType> spanningAlignments;

  /// List of indel sets, one per small alignment
  std::vector<CandidateSegmentSetType> smallSVSegments;

  /// Data specific to searching for a large insertion candidate
  std::vector<LargeInsertionInfo> largeInsertInfo;

  /// Remote reads retrieved to improve assembly and scoring for this locus
  RemoteReadCache remoteReads;

  /// If non-empty sv candidate set, which contig/alignment produced them?
  unsigned bestAlignmentIndex = 0;

  // expanded reference regions around the candidate SV breakend regions, for small events we use only bp1ref:
  reference_contig_segment bp1ref;
  reference_contig_segment bp2ref;

  std::vector<SVCandidate> svs;  ///< summarize candidate refined sv candidates

  /// If true, assembly was skipped for this case because of an overlapping assembly
  bool isOverlapSkip = false;
};
