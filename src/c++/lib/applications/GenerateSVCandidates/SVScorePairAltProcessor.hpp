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

#include "SVScorePairProcessor.hpp"
#include "manta/SVCandidateAssemblyData.hpp"
#include "manta/ShadowReadFinder.hpp"
#include "options/ReadScannerOptions.hpp"
#include "options/SVRefinerOptions.hpp"

struct ContigParams {
  ContigParams(const SVCandidateAssemblyData& assemblyData, const SVCandidate& sv);

  /// extended contig:
  const std::string& extSeq;

  /// where does the sv segment begin,end in reference coordinates?:
  known_pos_range2 segmentSpan;

  known_pos_range2 bpAOffset;
  known_pos_range2 bpBOffset;
};

/// estimate pair support for an sv candidate
/// restricted to simple indel style svs
struct SVScorePairAltProcessor : public SVScorePairProcessor {
  SVScorePairAltProcessor(
      const bam_header_info&         header,
      const ReadScannerOptions&      scanOpt,
      const SVRefinerOptions&        refineOpt,
      const std::vector<bool>&       initIsAlignmentTumor,
      const SVLocusScanner&          initReadScanner,
      const PairOptions&             initPairOpt,
      const SVCandidateAssemblyData& initAssemblyData,
      const SVCandidate&             initSv,
      const bool                     initIsBp1,
      SVEvidence&                    initEvidence)
    : SVScorePairProcessor(
          initIsAlignmentTumor, initReadScanner, initPairOpt, initSv, initIsBp1, initEvidence),
      _header(header),
      _assemblyData(initAssemblyData),
      _shadowAligner(refineOpt.spanningAlignScores),
      _shadow(
          scanOpt.minSingletonMapqCandidates,
          (!initIsBp1),  /// search for left-open shadows
          (initIsBp1)),  /// search for right-open shadows
      _contig(initAssemblyData, initSv)
  {
    checkInput(sv);
  }

  /// what to skip in addition to the core skip test?
  ///
  /// This override accounts for shadow and chimera re-maps for large insertions:
  ///
  bool isSkipRecord(const bam_record& bamRead) override
  {
    if (!isLargeInsertSV(sv)) return SVScorePairProcessor::isSkipRecord(bamRead);

    if (!bamRead.is_paired())
      return true;
    else if (bamRead.is_unmapped() && bamRead.is_mate_unmapped())
      return true;
    return false;
  }

  void processClearedRecord(
      const SVId& svId, const bam_record& bamRead, SVEvidenceWriterSampleData& svSupportFrags) override;

  /// TestSVScorerAltProcessor is a friend structure of SVScorePairAltProcessor. So that it can access private
  /// methods of SVScorePairAltProcessor.
  friend struct TestSVScorerAltProcessor;

private:
  static void checkInput(const SVCandidate& sv);

  /// \brief Realign one end of read pair
  ///
  /// \param[in] bamHeader This is only used for debug/exception messages.
  /// \param[in] fragmentQname QNAME of the anchor/floating read pair fragment. This is only used for
  /// debug/exception messages.
  /// \param[in] isLeftOfInsert If true, the anchor is on the left side of the insertion (and otherwise on the
  /// right)
  /// \param[in] floatRead The read to be realigned, already reverse complemented to the expected orientation,
  /// must not be empty
  /// \param[in] anchorTid The alignment contig id of the anchoring (ie. non-realigned) read of the pair. This
  /// is only used for debug/exception messages.
  /// \param[in] anchorPos The alignment position of the anchoring (ie. non-realigned) read of the pair
  /// \param[out] altTemplateSize The updated template size of the read pair following realignment.
  ///
  /// \return true for usable alignment
  bool realignPairedRead(
      const bam_header_info& bamHeader,
      const std::string&     fragmentQname,
      const bool             isLeftOfInsert,
      const std::string&     floatRead,
      const int              anchorTid,
      const pos_t            anchorPos,
      int&                   altTemplateSize);

  bool alignShadowRead(const bam_record& bamRead, int& altTemplateSize);

  /// test whether a frag reference span provides sufficient support for a breakpoint of this sv:
  bool testFragOverlap(const int fragBeginRefPos, const int fragEndRefPos) const;

  ///////////////////////
  const bam_header_info&         _header;
  const SVCandidateAssemblyData& _assemblyData;

  const GlobalAligner<int> _shadowAligner;
  ShadowReadFinder         _shadow;

  ContigParams _contig;
};
