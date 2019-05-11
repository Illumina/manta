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

#pragma once

#include "assembly/AssemblyReadInfo.hpp"
#include "blt_util/time_util.hpp"
#include "htsapi/bam_streamer.hpp"
#include "manta/ChromDepthFilterUtil.hpp"
#include "manta/SVCandidate.hpp"
#include "manta/SVCandidateAssemblyData.hpp"
#include "manta/SVLocusScanner.hpp"
#include "options/AlignmentFileOptions.hpp"
#include "options/IterativeAssemblerOptions.hpp"

#include <string>
#include <vector>

/// Assembles reads across suspected breakends of both 'spanning' (2 breakend regions) and
/// 'complex' (1 breakend region) low-resolution SVCandidates
///
/// In each case, the breakend regions are scanned for 'assembly-evidence' reads
/// consistent with the low-resolution SVCandidate hypothesis. These reads are assembled,
/// iterating over a range of word lengths until the first successful assembly. If unused
/// reads remain, the assembly is re-started using this subset.
///
struct SVCandidateAssembler {
  typedef IterativeAssemblerOptions AssemblerOptions;

  SVCandidateAssembler(
      const ReadScannerOptions&   scanOpt,
      const AssemblerOptions&     assembleOpt,
      const AlignmentFileOptions& alignFileOpt,
      const std::string&          referenceFilename,
      const std::string&          statsFilename,
      const std::string&          chromDepthFilename,
      const bam_header_info&      bamHeader,
      const AllSampleReadCounts&  counts,
      const bool                  isRNA,
      TimeTracker&                remoteReadRetrievalTime);

  /// Given a 'complex' SV candidate with 1 breakend region, assemble reads
  /// over the breakend region
  void assembleComplexSVCandidate(
      const SVBreakend&               bp,
      const reference_contig_segment& refSeq,
      const bool                      isSearchRemoteInsertionReads,
      RemoteReadCache&                remoteReads,
      Assembly&                       as) const;

  /// Given a 'spanning' SV candidate with 2 breakend regions, assemble reads
  /// over the junction of the 2 breakend regions
  void assembleSpanningSVCandidate(
      const SVBreakend&               bp1,
      const SVBreakend&               bp2,
      const bool                      isBp1Reversed,
      const bool                      isBp2Reversed,
      const reference_contig_segment& refSeq1,
      const reference_contig_segment& refSeq2,
      Assembly&                       as) const;

  const AssemblerOptions& getAssembleOpt() const { return _assembleOpt; }

  typedef std::map<std::string, unsigned> ReadIndexType;

private:
  typedef std::shared_ptr<bam_streamer> streamPtr;

  /// Collect reads crossing an SV breakpoint and add them to 'reads'
  ///
  /// \param[in] isReversed if true revcomp all reads on input
  ///
  /// \param[in] refSeq this is used to find reads which have poorly aligned ends, such reads are added to the
  /// breakend assembly pool
  ///
  /// \param[in] isSearchRemoteInsertionReads if true search the remote end of chimeric pairs for MAPQ0
  /// insertion support
  ///
  /// \param[out] remoteReadsCache stores any discovered remote reads so that these can be reused during
  /// scoring
  ///
  /// \param[out] reads collected breakend assembly candidate reads
  void getBreakendReads(
      const SVBreakend&               bp,
      const bool                      isReversed,
      const reference_contig_segment& refSeq,
      const bool                      isSearchRemoteInsertionReads,
      RemoteReadCache&                remoteReadsCache,
      ReadIndexType&                  readIndex,
      AssemblyReadInput&              reads) const;

  const ReadScannerOptions _scanOpt;
  const AssemblerOptions   _assembleOpt;
  const std::vector<bool>  _isAlignmentTumor;

  /// This depth filters tracks the maximum depth for procssing and assembling a locus
  const ChromDepthFilterUtil _dFilter;

  /// This depth filter tracks the maximum depth of a source locus which qualifies for remote read recovery
  /// used to assemble large insertions
  const ChromDepthFilterUtil _dFilterLocalDepthForRemoteReadRetrieval;

  // contains functions to detect/classify anomalous reads
  SVLocusScanner         _readScanner;
  std::vector<streamPtr> _bamStreams;
  TimeTracker&           _remoteReadRetrievalTime;

  /// In each sample, store the background rate of reads which would qualify for remove recovery at an
  /// insertion locus
  std::vector<double> _sampleRemoteRecoveryCandidateRate;
};
