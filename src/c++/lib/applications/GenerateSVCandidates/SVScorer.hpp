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

#pragma once

#include <string>
#include <vector>

#include "GSCOptions.hpp"
#include "JunctionCallInfo.hpp"
#include "SVEvidence.hpp"
#include "SVScorePairProcessor.hpp"
#include "SplitReadAlignment.hpp"
#include "assembly/AssembledContig.hpp"
#include "blt_util/qscore_snp.hpp"
#include "htsapi/bam_header_info.hpp"
#include "htsapi/bam_streamer.hpp"
#include "manta/ChromDepthFilterUtil.hpp"
#include "manta/SVCandidateAssemblyData.hpp"
#include "manta/SVCandidateSetData.hpp"
#include "manta/SVLocusScanner.hpp"
#include "manta/SVModelScoreInfo.hpp"
#include "manta/SVMultiJunctionCandidate.hpp"
#include "manta/SVScoreInfoSomatic.hpp"

struct CallOptionsSharedDeriv {
  explicit CallOptionsSharedDeriv(const CallOptionsShared& opt) : refQ(opt.snpPrior), altQ(0) {}

  const qscore_snp refQ;
  const qscore_snp altQ;
};

struct CallOptionsDiploidDeriv : private boost::noncopyable {
  explicit CallOptionsDiploidDeriv(const CallOptionsDiploid& opt)
  {
    using namespace DIPLOID_GT;

    assert(opt.indelPrior < 0.5);

    prior[HET] = opt.indelPrior;
    prior[HOM] = opt.indelPrior / 2;
    prior[REF] = 1 - prior[HET] - prior[HOM];

    for (unsigned i(0); i < SIZE; ++i) {
      logPrior[i] = std::log(prior[i]);
    }
  }

  std::array<float, DIPLOID_GT::SIZE> prior;
  std::array<float, DIPLOID_GT::SIZE> logPrior;
};

struct CallOptionsSomaticDeriv : private boost::noncopyable {
  explicit CallOptionsSomaticDeriv(const CallOptionsSomatic& opt)
  {
    using namespace SOMATIC_GT;

    assert(opt.germlineSVPrior < 0.5);

    prior[SOM]   = opt.somaticSVPrior;
    prior[NOISE] = opt.largeNoiseSVPrior;

    prior[HET] = opt.germlineSVPrior;
    prior[HOM] = opt.germlineSVPrior / 2;

    // this assumes all states independent, and somatic and noise only occur on germline ref GT background:
    const float nonref(prior[SOM] + prior[NOISE] + prior[HET] + prior[HOM]);
    assert(nonref >= 0 && nonref <= 1);
    prior[REF] = 1 - nonref;

    for (unsigned i(0); i < SIZE; ++i) {
      _logPrior[i] = std::log(prior[i]);
    }

    smallNoisePrior    = opt.smallNoiseSVPrior;
    largeNoisePrior    = opt.largeNoiseSVPrior;
    logSmallNoisePrior = std::log(opt.smallNoiseSVPrior);
    logLargeNoisePrior = std::log(opt.largeNoiseSVPrior);
  }

  float logPrior(const unsigned gt, const float largeNoiseWeight) const
  {
    assert(largeNoiseWeight >= 0. && largeNoiseWeight <= 1.);

    if (gt != SOMATIC_GT::NOISE) return _logPrior[gt];

    if (largeNoiseWeight <= 0.) return logSmallNoisePrior;
    if (largeNoiseWeight >= 1.) return logLargeNoisePrior;

    return std::log((1 - largeNoiseWeight) * smallNoisePrior + largeNoiseWeight * largeNoisePrior);
  }

private:
  std::array<float, SOMATIC_GT::SIZE> prior;
  std::array<float, SOMATIC_GT::SIZE> _logPrior;

  float smallNoisePrior;
  float largeNoisePrior;
  float logSmallNoisePrior;
  float logLargeNoisePrior;
};

/// Implements SV scoring/genotyping process
///
struct SVScorer {
  SVScorer(const GSCOptions& opt, const SVLocusScanner& readScanner, const bam_header_info& header);

  /// Gather supporting read evidence and generate:
  /// 1. diploid quality score and genotype for SV candidate
  /// 2. somatic quality score
  ///
  /// \param mjModelScoreInfo standard (junction-independent) scoring information for each junction
  ///
  /// \param mjJointModelScoreInfo  joint junction scoring information used for cases where a multi-junction
  /// event is detected
  void scoreSV(
      const SVCandidateSetData&                   svData,
      const std::vector<SVCandidateAssemblyData>& mjAssemblyData,
      const SVMultiJunctionCandidate&             mjSV,
      const std::vector<SVId>&                    mjSVId,
      const std::vector<bool>&                    isJunctionFiltered,
      const bool                                  isSomatic,
      const bool                                  isTumorOnly,
      std::vector<SVModelScoreInfo>&              mjModelScoreInfo,
      SVModelScoreInfo&                           mjJointModelScoreInfo,
      bool&                                       isMJEvent,
      SVEvidenceWriterData&                       svEvidenceWriterData);

  typedef std::shared_ptr<SVScorePairProcessor> pairProcPtr;
  typedef std::shared_ptr<bam_streamer>         streamPtr;

#if 0
    unsigned
    sampleCount() const
    {
        return _sampleCount;
    }

    const std::vector<std::string>&
    sampleNames() const
    {
        return _sampleNames;
    }
#endif

  /// TestSVScorer is a friend structure of SV scorer. So that it can access private
  /// methods of SVScorer. As SVScorer has many private methods, for unit test writing
  /// this friend structure has been created.
  friend struct TestSVScorer;

private:
  void processExistingAltPairInfo(
      const PairOptions&        pairOpt,
      const SVCandidateSetData& svData,
      const SVCandidate&        sv,
      const SVId&               svId,
      SVEvidence&               evidence,
      SVEvidenceWriterData&     svSupports);

  /// Estimate pair support for an sv candidate
  ///
  /// This is restricted to simple indel-style SVs
  void getSVAltPairSupport(
      const PairOptions&             pairOpt,
      const SVCandidateAssemblyData& assemblyData,
      const SVCandidate&             sv,
      SVEvidence&                    evidence,
      std::vector<pairProcPtr>&      pairProcList);

  /// Find spanning read support for the reference allele for sv candidate
  void getSVRefPairSupport(
      const PairOptions&        pairOpt,
      const SVCandidate&        sv,
      SVEvidence&               evidence,
      std::vector<pairProcPtr>& pairProcList);

  /// Find paired read support for ref and alt alleles
  void getSVPairSupport(
      const SVCandidateSetData&      svData,
      const SVCandidateAssemblyData& assemblyData,
      const SVCandidate&             sv,
      const SVId&                    svId,
      SVEvidence&                    evidence,
      SVEvidenceWriterData&          svSupports);

  /// \brief Find split read support for ref and alt alleles
  ///
  /// \param assemblyData Assembly results for the given candidate SV
  ///
  /// \param sv Candidate SV for which split read support is being enumerated
  ///
  /// \param svId Structural variant ID tag used to annotate BAM output used for debugging
  ///
  /// \param ssInfo Model-agnostic summary information for the candidate SV
  ///
  /// \param evidence The structure used to return all split read information for the variant, used for
  /// generating scores under specific variant models.
  ///
  /// \param svSupports Structure into which supporting evidence reads can be added for BAM output used in
  /// debugging
  void getSVSplitReadSupport(
      const SVCandidateAssemblyData& assemblyData,
      const SVCandidate&             sv,
      const SVId&                    svId,
      SVScoreInfo&                   ssInfo,
      SVEvidence&                    evidence,
      SVEvidenceWriterData&          svSupports);

  /// Determine maximum depth and MQ0 frac in region around breakend of normal sample
  void getBreakendMaxMappedDepthAndMQ0(
      const bool        isTumorOnly,
      const bool        isMaxDepth,
      const double      cutoffDepth,
      const SVBreakend& bp,
      unsigned&         maxDepth,
      float&            MQ0Frac);

  /// Apply all scoring models relevant to this event:
  ///
  /// \param junctionData one element describing each junction of an event, for normal (single-junction)
  /// candidates, the vector size should be one
  void computeAllScoreModels(
      const bool                           isSomatic,
      const bool                           isTumorOnly,
      const std::vector<JunctionCallInfo>& junctionData,
      SVModelScoreInfo&                    modelScoreInfo);

  /// Accumulate (model-agnostic) evidence for the SV alt/ref alleles
  ///
  /// This step gathers information (such as counts and allele likelihoods) to be used downstream by more
  /// specific (germline, somatic, etc..) scoring models.
  ///
  void getSVSupportingEvidence(
      const SVCandidateSetData&      svData,
      const SVCandidateAssemblyData& assemblyData,
      const bool                     isTumorOnly,
      const SVCandidate&             sv,
      const SVId&                    svId,
      SVScoreInfo&                   ssInfo,
      SVEvidence&                    evidence,
      SVEvidenceWriterData&          svSupports);

  const std::vector<bool>       _isAlignmentTumor;
  const bool                    _isRNA;
  const CallOptionsShared       _callOpt;
  const CallOptionsSharedDeriv  _callDopt;
  const CallOptionsDiploid      _diploidOpt;
  const CallOptionsDiploidDeriv _diploidDopt;
  const ReadScannerOptions      _scanOpt;
  const SVRefinerOptions        _refineOpt;
  const CallOptionsSomatic      _somaticOpt;
  const CallOptionsSomaticDeriv _somaticDopt;
  const CallOptionsTumor        _tumorOpt;
  const ChromDepthFilterUtil    _dFilterDiploid;
  const ChromDepthFilterUtil    _dFilterSomatic;
  const ChromDepthFilterUtil    _dFilterTumor;
  const SVLocusScanner&         _readScanner;

  const bam_header_info& _header;

  std::vector<streamPtr> _bamStreams;

  unsigned _sampleCount;
  unsigned _diploidSampleCount;
};
