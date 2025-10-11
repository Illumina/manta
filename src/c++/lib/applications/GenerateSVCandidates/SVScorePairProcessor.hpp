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

#include "SVEvidence.hpp"
#include "SVEvidenceWriter.hpp"
#include "SVScorerPairOptions.hpp"
#include "SVScorerShared.hpp"
#include "blt_util/SizeDistribution.hpp"
#include "manta/BamRegionProcessor.hpp"
#include "manta/JunctionIdGenerator.hpp"
#include "manta/ReadFilter.hpp"
#include "manta/SVCandidate.hpp"
#include "manta/SVLocusScanner.hpp"

struct SVScorePairInitParams {
  SVScorePairInitParams(const SVLocusScanner& readScanner, const SVCandidate& sv, const bool isBp1);

  pos_t centerPosA;
  pos_t centerPosB;
  pos_t centerPos;

  // total impact of the alt allele on template size:
  pos_t    altShift;
  unsigned minMapQ;
  unsigned minTier2MapQ;  // a second, lower mapq threshold used to disprove a somatic allele during
                          // tumor/normal calling
};

struct SVScorePairBamParams {
  bool                    isSet         = false;
  unsigned                bamIndex      = 0;
  bool                    isTumor       = false;
  pos_t                   minFrag       = 0;
  pos_t                   maxFrag       = 0;
  const SizeDistribution* fragDistroPtr = nullptr;
  GenomeInterval          interval;
};

struct SVScorePairProcessor : public BamRegionProcessor {
  SVScorePairProcessor(
      const std::vector<bool>& initIsAlignmentTumor,
      const SVLocusScanner&    initReadScanner,
      const PairOptions&       initPairOpt,
      const SVCandidate&       initSv,
      const bool               initIsBp1,
      SVEvidence&              initEvidence)
    : isAlignmentTumor(initIsAlignmentTumor),
      readScanner(initReadScanner),
      pairOpt(initPairOpt),
      sv(initSv),
      isBp1(initIsBp1),
      evidence(initEvidence),
      svParams(readScanner, sv, isBp1),
      bamParams()
  {
  }

  const GenomeInterval& nextBamIndex(const unsigned bamIndex);

  // alternate interface
  static bool isSkipRecordCore(const bam_record& bamRead)
  {
    return (isReadFilteredCore(bamRead) || bamRead.isNonStrictSupplement());
  }

  /// what to skip in addition to the core skip test?
  virtual bool isSkipRecord(const bam_record& bamRead)
  {
    if (bamRead.is_unmapped() || (bamRead.is_paired() && bamRead.is_mate_unmapped()))
      return true;
    else if (!is_innie_pair(bamRead))
      return true;
    return false;
  }

  // process a record for which isSkipRecord() == false
  virtual void processClearedRecord(
      const SVId& svId, const bam_record& bamRead, SVEvidenceWriterSampleData& svSupportFrags) = 0;

  static bool isLargeInsertSV(const SVCandidate& sv) { return (sv.insertSeq.size() >= 100); }

protected:
  void setAlleleFrag(
      const SizeDistribution&           fragDistro,
      const int                         size,
      SVFragmentEvidenceAlleleBreakend& bp,
      const bool /*isPdf*/ = false)
  {
    float fragProb(0);
#if 0
        if (isPdf)
        {
            fragProb = fragDistro.pdf(size);
        }
        else
#endif
    {
      fragProb = fragDistro.cdf(size);
      fragProb = std::min(fragProb, (1 - fragProb));
    }
    if (pairOpt.RNA) {
      fragProb = std::max(fragProb, pairOpt.minFragProb);
    }
#ifdef DEBUG_MEGAPAIR
    log_os << __FUNCTION__ << ": fraglen,prob " << size << " " << fragProb << "\n";
#endif
    bp.isFragmentSupport = true;
    bp.fragLengthProb    = fragProb;
  }

  const std::vector<bool> isAlignmentTumor;
  const SVLocusScanner&   readScanner;
  const PairOptions&      pairOpt;
  const SVCandidate&      sv;
  const bool              isBp1;
  SVEvidence&             evidence;

  const SVScorePairInitParams svParams;
  SVScorePairBamParams        bamParams;
};
