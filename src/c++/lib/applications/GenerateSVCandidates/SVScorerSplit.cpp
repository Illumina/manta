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
/// \author Chris Saunders and Xiaoyu Chen
///
/// This file contains a subset of the implementation for SVScorer.hpp focusing on split read handling
///

#include "SVScorer.hpp"

#include <iostream>

#include "boost/scoped_array.hpp"

#include "blt_util/log.hpp"
#include "blt_util/seq_util.hpp"
#include "htsapi/SimpleAlignment_bam_util.hpp"
#include "manta/ReadFilter.hpp"
#include "manta/ShadowReadFinder.hpp"

//#define DEBUG_SVS

//#define DEBUG_SUPPORT

static void incrementAlleleEvidence(
    const SRAlignmentInfo&                   bp1SR,
    const SRAlignmentInfo&                   bp2SR,
    const unsigned                           readMapQ,
    const float                              evidence,
    SVSampleAlleleInfo&                      allele,
    SVFragmentEvidenceAlleleBreakendPerRead& bp1Support,
    SVFragmentEvidenceAlleleBreakendPerRead& bp2Support)
{
  if (bp1SR.isEvidence) {
    bp1Support.isSplitSupport = true;
    bp1Support.splitEvidence  = bp1SR.evidence;
  }

  if (bp1SR.isTier2Evidence) {
    bp1Support.isTier2SplitSupport = true;
  }

  if (bp2SR.isEvidence) {
    bp2Support.isSplitSupport = true;
    bp2Support.splitEvidence  = bp2SR.evidence;
  }

  if (bp2SR.isTier2Evidence) {
    bp2Support.isTier2SplitSupport = true;
  }

  if (bp1SR.isEvidence || bp2SR.isEvidence) {
    allele.splitReadCount++;
    allele.splitReadEvidence += evidence;
    allele.splitReadMapQ += readMapQ * readMapQ;

#ifdef DEBUG_SVS
    log_os << __FUNCTION__ << "\n";
    log_os << "bp1\n";
    log_os << bp1SR;
    log_os << "bp2\n";
    log_os << bp2SR;
    log_os << "evidence = " << evidence << "\n";
    log_os << "accumulated evidence = " << allele.splitReadEvidence << "\n";
    log_os << "split read Count = " << allele.splitReadCount << "\n\n";
#endif
  }
}

static void incrementSplitReadEvidence(
    const SRAlignmentInfo&                   refBp1SR,
    const SRAlignmentInfo&                   refBp2SR,
    const SRAlignmentInfo&                   altBp1SR,
    const SRAlignmentInfo&                   altBp2SR,
    const unsigned                           readMapQ,
    const bool                               isRNA,
    SVSampleAlleleInfo&                      refAllele,
    SVSampleAlleleInfo&                      altAllele,
    SVFragmentEvidenceAlleleBreakendPerRead& refBp1Support,
    SVFragmentEvidenceAlleleBreakendPerRead& refBp2Support,
    SVFragmentEvidenceAlleleBreakendPerRead& altBp1Support,
    SVFragmentEvidenceAlleleBreakendPerRead& altBp2Support)
{
  refBp1Support.splitLnLhood = refBp1SR.alignLnLhood;
  refBp2Support.splitLnLhood = refBp2SR.alignLnLhood;
  altBp1Support.splitLnLhood = altBp1SR.alignLnLhood;
  altBp2Support.splitLnLhood = altBp2SR.alignLnLhood;

  const float refAlignLlh(std::max(refBp1SR.alignLnLhood, refBp2SR.alignLnLhood));
  const float altAlignLlh(std::max(altBp1SR.alignLnLhood, altBp2SR.alignLnLhood));

  // For DNA, a split read is considered supporting the breakend
  // if the read's alignment to the allele with higher likelihood meets the evidence requirements
  if (isRNA || (refAlignLlh > altAlignLlh)) {
    float       refBp1Evidence = refBp1SR.isEvidence ? refBp1SR.evidence : 0;
    float       refBp2Evidence = refBp2SR.isEvidence ? refBp2SR.evidence : 0;
    const float refEvidence(std::max(refBp1Evidence, refBp2Evidence));

    incrementAlleleEvidence(
        refBp1SR, refBp2SR, readMapQ, refEvidence, refAllele, refBp1Support, refBp2Support);
  }

  if (isRNA || (altAlignLlh > refAlignLlh)) {
    float       altBp1Evidence = altBp1SR.isEvidence ? altBp1SR.evidence : 0;
    float       altBp2Evidence = altBp2SR.isEvidence ? altBp2SR.evidence : 0;
    const float altEvidence(std::max(altBp1Evidence, altBp2Evidence));

    incrementAlleleEvidence(
        altBp1SR, altBp2SR, readMapQ, altEvidence, altAllele, altBp1Support, altBp2Support);
  }
}

/// Score a single split read for one breakend of a candidate SV in one sample
static void getReadSplitScore(
    const CallOptionsSharedDeriv&   dopt,
    const unsigned                  flankScoreSize,
    const SVId&                     svId,
    const SVBreakend&               bp,
    const SVAlignmentInfo&          svAlignInfo,
    const reference_contig_segment& bpRef,
    const bool                      isBP1,
    const unsigned                  minMapQ,
    const unsigned                  minTier2MapQ,
    const bool                      isRNA,
    const bool                      isShadow,
    const bool                      isReversedShadow,
    const bam_record&               bamRead,
    SVEvidence::evidenceTrack_t&    sampleEvidence,
    SVSampleInfo&                   sample,
    SVEvidenceWriterSampleData&     svSupportFrags)
{
#ifdef DEBUG_SVS
  log_os << __FUNCTION__ << " split scoring read: " << bamRead << "\n";
#endif

  SVFragmentEvidence& fragment(sampleEvidence[bamRead.qname()]);

  const bool                               isRead1(bamRead.is_first());
  SVFragmentEvidenceAlleleBreakendPerRead& altBp1ReadSupport(fragment.alt.bp1.getRead(isRead1));

  // In this function we evaluate the hypothesis of both breakends at the same time, the only difference bp1
  // vs bp2 makes is where in the bam we look for reads, therefore if we see split evaluation for bp1 or bp2,
  // we can skip this read:
  if (altBp1ReadSupport.isSplitEvaluated) return;

  SVFragmentEvidenceAlleleBreakendPerRead& refBp1ReadSupport(fragment.ref.bp1.getRead(isRead1));
  SVFragmentEvidenceAlleleBreakendPerRead& altBp2ReadSupport(fragment.alt.bp2.getRead(isRead1));
  SVFragmentEvidenceAlleleBreakendPerRead& refBp2ReadSupport(fragment.ref.bp2.getRead(isRead1));
  altBp1ReadSupport.isSplitEvaluated = true;
  refBp1ReadSupport.isSplitEvaluated = true;
  altBp2ReadSupport.isSplitEvaluated = true;
  refBp2ReadSupport.isSplitEvaluated = true;

  std::string    readSeq = bamRead.get_bam_read().get_string();
  const uint8_t* qual(bamRead.qual());

  // Reverse read and qual sequence if this is a reversed shadow read:
  boost::scoped_array<uint8_t> qualcpy;
  if (isShadow && isReversedShadow) {
    reverseCompStr(readSeq);

    qualcpy.reset(new uint8_t[readSeq.size()]);
    std::reverse_copy(qual, qual + readSeq.size(), qualcpy.get());
    qual = qualcpy.get();
  }

  SVFragmentEvidenceRead& evidenceRead(fragment.getRead(isRead1));
  setReadEvidence(minMapQ, minTier2MapQ, bamRead, isShadow, evidenceRead);

  // align the read to the alt allele contig
  SRAlignmentInfo altBp1SR;
  SRAlignmentInfo altBp2SR;
  splitReadAligner(
      flankScoreSize,
      readSeq,
      dopt.altQ,
      qual,
      svAlignInfo.bp1ContigSeq(),
      svAlignInfo.bp1ContigOffset,
      altBp1SR);
  splitReadAligner(
      flankScoreSize,
      readSeq,
      dopt.altQ,
      qual,
      svAlignInfo.bp2ContigSeq(),
      svAlignInfo.bp2ContigOffset,
      altBp2SR);

  // align the read to reference regions
  SRAlignmentInfo refBp1SR;
  SRAlignmentInfo refBp2SR;
  if (!isRNA) {
    splitReadAligner(
        flankScoreSize,
        readSeq,
        dopt.refQ,
        qual,
        svAlignInfo.bp1ReferenceSeq(),
        svAlignInfo.bp1RefOffset,
        refBp1SR);
    splitReadAligner(
        flankScoreSize,
        readSeq,
        dopt.refQ,
        qual,
        svAlignInfo.bp2ReferenceSeq(),
        svAlignInfo.bp2RefOffset,
        refBp2SR);
  } else {
    if (isBP1)
      getRefAlignment(bamRead, bpRef, bp.interval.range, dopt.refQ, refBp1SR);
    else
      getRefAlignment(bamRead, bpRef, bp.interval.range, dopt.refQ, refBp2SR);
  }
#ifdef DEBUG_SVS
  log_os << "\t reference align bp1: " << refBp1SR << "\n";
  log_os << "\t reference align bp2: " << refBp2SR << "\n";
#endif

  const unsigned readMapQ = bamRead.map_qual();
  // scoring
  incrementSplitReadEvidence(
      refBp1SR,
      refBp2SR,
      altBp1SR,
      altBp2SR,
      readMapQ,
      isRNA,
      sample.ref,
      sample.alt,
      refBp1ReadSupport,
      refBp2ReadSupport,
      altBp1ReadSupport,
      altBp2ReadSupport);

  if (fragment.isAltSplitReadSupport(bamRead.is_first())) {
    SVEvidenceWriterReadPair& supportFrag(svSupportFrags.getSupportFragment(bamRead));
    supportFrag.addSplitSupport(bamRead.is_first(), svId.localId);

#ifdef DEBUG_SUPPORT
    log_os << __FUNCTION__ << "  Adding a split read that supports the alt allele: " << bamRead.qname();
    if (bamRead.is_first())
      log_os << "\tR1";
    else
      log_os << "\tR2";
    log_os << "\n" << supportFrag;
#endif
  }
}

/// Score split reads for one breakend of a candidate SV in one sample
///
/// \param flankScoreSize Number of bases from the end of the microhomology range used for the split read
/// score.
///
/// \param svId Structural variant ID tag used to annotate BAM output used for debugging
///
/// \param bp The breakend targeted for split read evaluation
///
/// \param svAlignInfo Details how the breakend maps to sv contig and reference
///
static void scoreSplitReads(
    const CallOptionsSharedDeriv&   dopt,
    const unsigned                  flankScoreSize,
    const SVId&                     svId,
    const SVBreakend&               bp,
    const SVAlignmentInfo&          svAlignInfo,
    const reference_contig_segment& bpRef,
    const bool                      isBP1,
    const unsigned                  minMapQ,
    const unsigned                  minTier2MapQ,
    const int                       bamShadowSearchDistance,
    const unsigned                  shadowMinMapq,
    const bool                      isRNA,
    SVEvidence::evidenceTrack_t&    sampleEvidence,
    bam_streamer&                   readStream,
    SVSampleInfo&                   sample,
    SVEvidenceWriterSampleData&     svSupportFrags)
{
  // Window to look for alignments that may (if unclipped) overlap the breakpoint:
  static const int extendedSearchRange(200);

  // Extract reads overlapping the break point.
  //
  // We are not looking for remote reads, (semialigned-) reads mapping near this breakpoint, but not across it
  // or any other kind of additional reads used for assembly.

  readStream.resetRegion(
      bp.interval.tid,
      std::max(0, bp.interval.range.begin_pos() - extendedSearchRange),
      bp.interval.range.end_pos() + extendedSearchRange);

  while (readStream.next()) {
    const bam_record& bamRead(*(readStream.get_record_ptr()));

    if (isReadUnmappedOrFilteredCore(bamRead)) continue;

    // TODO: remove this filter?
    // The supplemental alignment is likely to be hard-clipped
    if (bamRead.isNonStrictSupplement()) continue;

    // Skip reads that do not overlap the entire homology range of this breakpoint.
    const known_pos_range2 bamRange(matchifyEdgeSoftClipRefRange(getAlignment(bamRead)));
    if (!bamRange.is_range_intersect(bp.interval.range)) continue;

    static const bool isShadow(false);
    static const bool isReversedShadow(false);

    try {
      getReadSplitScore(
          dopt,
          flankScoreSize,
          svId,
          bp,
          svAlignInfo,
          bpRef,
          isBP1,
          minMapQ,
          minTier2MapQ,
          isRNA,
          isShadow,
          isReversedShadow,
          bamRead,
          sampleEvidence,
          sample,
          svSupportFrags);
    } catch (...) {
      log_os << "ERROR: Exception caught in getReadSplitScore() while scoring read: " << bamRead << "\n";
      throw;
    }
  }

  static const bool isIncludeShadowReads(false);

  // search for appropriate shadow reads to add to the split read pool
  //
  if (isIncludeShadowReads) {
    // depending on breakend type we may only be looking for candidates in one direction:
    bool             isSearchForLeftOpen(true);
    bool             isSearchForRightOpen(true);
    known_pos_range2 shadowRange;
    if (bp.state == SVBreakendState::RIGHT_OPEN) {
      isSearchForLeftOpen = false;

      shadowRange.set_begin_pos(std::max(0, bp.interval.range.begin_pos() - bamShadowSearchDistance));
      shadowRange.set_end_pos(bp.interval.range.begin_pos());
    } else if (bp.state == SVBreakendState::LEFT_OPEN) {
      isSearchForRightOpen = false;

      shadowRange.set_begin_pos(bp.interval.range.end_pos());
      shadowRange.set_end_pos(bp.interval.range.end_pos() + bamShadowSearchDistance);
    } else {
      assert(false && "Invalid bp state");
    }

    readStream.resetRegion(bp.interval.tid, shadowRange.begin_pos(), shadowRange.end_pos());

    ShadowReadFinder shadow(shadowMinMapq, isSearchForLeftOpen, isSearchForRightOpen);

    while (readStream.next()) {
      const bam_record& bamRead(*(readStream.get_record_ptr()));

      if (isReadFilteredCore(bamRead)) continue;
      if (!shadow.check(bamRead)) continue;

      static const bool isShadow(true);
      const bool        isReversedShadow(bamRead.is_mate_fwd_strand());

      //const uint8_t mapq(shadow.getMateMapq());
      getReadSplitScore(
          dopt,
          flankScoreSize,
          svId,
          bp,
          svAlignInfo,
          bpRef,
          isBP1,
          minMapQ,
          minTier2MapQ,
          isRNA,
          isShadow,
          isReversedShadow,
          bamRead,
          sampleEvidence,
          sample,
          svSupportFrags);
    }
  }
}

/// return rms given sum of squares
static float finishRms(const float sumSqr, const unsigned count)
{
  if (count == 0) return 0.;
  return std::sqrt(sumSqr / static_cast<float>(count));
}

static void finishRms(SVSampleAlleleInfo& sai)
{
  sai.splitReadMapQ = finishRms(sai.splitReadMapQ, sai.splitReadCount);
}

/// make final split read computations after bam scanning is finished:
static void finishSampleSRData(SVSampleInfo& sample)
{
  // finish rms mapq:
  finishRms(sample.alt);
  finishRms(sample.ref);
}

void SVScorer::getSVSplitReadSupport(
    const SVCandidateAssemblyData& assemblyData,
    const SVCandidate&             sv,
    const SVId&                    svId,
    SVScoreInfo&                   baseInfo,
    SVEvidence&                    evidence,
    SVEvidenceWriterData&          svSupports)
{
  // Apply the split-read scoring only when:
  // 1. the SV is precise, i.e. has successfully aligned contigs;
  // 2. the values of max depth are reasonable (otherwise, the read map may blow out). Note the max depth
  // filter is run externally
  if (sv.isImprecise()) return;

  // Get data on standard read pairs crossing the two breakends,

  // extract SV alignment info for split read evidence
  const SVAlignmentInfo SVAlignInfo(sv, assemblyData);

  // Number of bases from the end of the microhomology range used for the split read score:
  static const unsigned flankScoreSize(50);

  // only consider a split alignment with sufficient flanking sequence:
  static const unsigned minBreakpointDistanceFromEdge(100);
  if (!SVAlignInfo.isMinBpEdge(minBreakpointDistanceFromEdge)) return;

#ifdef DEBUG_SVS
  log_os << __FUNCTION__ << " sv: " << sv << '\n';
  log_os << __FUNCTION__ << " SVAlignInfo: " << SVAlignInfo << '\n';
#endif

  const unsigned minMapQ(_readScanner.getMinMapQ());
  const unsigned minTier2MapQ(_readScanner.getMinTier2MapQ());

  const unsigned bamCount(_bamStreams.size());
  for (unsigned bamIndex(0); bamIndex < bamCount; ++bamIndex) {
    SVSampleInfo& sample(baseInfo.samples[bamIndex]);
    bam_streamer& bamStream(*_bamStreams[bamIndex]);

    SVEvidence::evidenceTrack_t& sampleEvidence(evidence.getSampleEvidence(bamIndex));
    SVEvidenceWriterSampleData&  svSupportFrags(svSupports.getSampleData(bamIndex));

    const int bamShadowSearchDistance(_readScanner.getShadowSearchDistance(bamIndex));

    // scoring split reads overlapping bp1
#ifdef DEBUG_SVS
    log_os << __FUNCTION__ << " scoring BP1 " << sv.bp1.interval << "\n";
#endif
    scoreSplitReads(
        _callDopt,
        flankScoreSize,
        svId,
        sv.bp1,
        SVAlignInfo,
        assemblyData.bp1ref,
        true,
        minMapQ,
        minTier2MapQ,
        bamShadowSearchDistance,
        _scanOpt.minSingletonMapqCandidates,
        _isRNA,
        sampleEvidence,
        bamStream,
        sample,
        svSupportFrags);
    // scoring split reads overlapping bp2
#ifdef DEBUG_SVS
    log_os << __FUNCTION__ << " scoring BP2 " << sv.bp2.interval << "\n";
#endif
    scoreSplitReads(
        _callDopt,
        flankScoreSize,
        svId,
        sv.bp2,
        SVAlignInfo,
        assemblyData.bp2ref,
        false,
        minMapQ,
        minTier2MapQ,
        bamShadowSearchDistance,
        _scanOpt.minSingletonMapqCandidates,
        _isRNA,
        sampleEvidence,
        bamStream,
        sample,
        svSupportFrags);

    finishSampleSRData(sample);
  }

#ifdef DEBUG_SVS
  for (unsigned bamIndex(0); bamIndex < bamCount; ++bamIndex) {
    log_os << "bam index: " << bamIndex << "\n";
    const SVSampleInfo& sample(baseInfo.samples[bamIndex]);
    log_os << "Alt contig SR count: " << sample.alt.splitReadCount << "\n";
    log_os << "Alt contig SR evidence: " << sample.alt.splitReadEvidence << "\n";
    log_os << "Alt contig SR_mapQ: " << sample.alt.splitReadMapQ << "\n";

    log_os << "Ref contig SR count: " << sample.ref.splitReadCount << "\n";
    log_os << "Ref contig SR evidence: " << sample.ref.splitReadEvidence << "\n";
    log_os << "Ref contig SR_mapQ: " << sample.ref.splitReadMapQ << "\n";
  }
#endif
}
