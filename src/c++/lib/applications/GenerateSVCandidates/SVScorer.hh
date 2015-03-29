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
/// \author Chris Saunders and Xiaoyu Chen
///

#pragma once

#include "GSCOptions.hh"
#include "JunctionCallInfo.hh"
#include "SplitReadAlignment.hh"
#include "SVEvidence.hh"
#include "SVScorerPairOptions.hh"
#include "SVScorePairProcessor.hh"

#include "assembly/AssembledContig.hh"
#include "blt_util/qscore_snp.hh"
#include "htsapi/bam_streamer.hh"
#include "htsapi/bam_header_info.hh"
#include "manta/ChromDepthFilterUtil.hh"
#include "manta/SVCandidateAssemblyData.hh"
#include "manta/SVCandidateSetData.hh"
#include "manta/SVLocusScanner.hh"
#include "manta/SVModelScoreInfo.hh"
#include "manta/SVMultiJunctionCandidate.hh"
#include "manta/SVScoreInfoSomatic.hh"

#include <vector>
#include <string>


struct CallOptionsSharedDeriv
{
    CallOptionsSharedDeriv(
        const CallOptionsShared& opt) :
        refQ(opt.snpPrior),
        altQ(0)
    {}

    const qscore_snp refQ;
    const qscore_snp altQ;
};


struct CallOptionsDiploidDeriv : private boost::noncopyable
{
    CallOptionsDiploidDeriv(
        const CallOptionsDiploid& opt)
    {
        using namespace DIPLOID_GT;

        assert(opt.indelPrior < 0.5);

        prior[HET] = opt.indelPrior;
        prior[HOM] = opt.indelPrior/2;
        prior[REF] = 1 - prior[HET] - prior[HOM];

        for (unsigned i(0); i<SIZE; ++i)
        {
            logPrior[i] = std::log(prior[i]);
        }
    }

    std::array<float,DIPLOID_GT::SIZE> prior;
    std::array<float,DIPLOID_GT::SIZE> logPrior;
};


struct CallOptionsSomaticDeriv : private boost::noncopyable
{
    CallOptionsSomaticDeriv(
        const CallOptionsSomatic& opt)
    {
        using namespace SOMATIC_GT;

        assert(opt.germlineSVPrior < 0.5);

        prior[SOM] = opt.somaticSVPrior;
        prior[NOISE] = opt.largeNoiseSVPrior;

        prior[HET] = opt.germlineSVPrior;
        prior[HOM] = opt.germlineSVPrior/2;

        // this assumes all states independent, and somatic and noise only occur on germline ref GT background:
        const float nonref(prior[SOM]+prior[NOISE]+prior[HET]+prior[HOM]);
        assert(nonref>=0 && nonref<=1);
        prior[REF] = 1 - nonref;

        for (unsigned i(0); i<SIZE; ++i)
        {
            _logPrior[i] = std::log(prior[i]);
        }

        smallNoisePrior = opt.smallNoiseSVPrior;
        largeNoisePrior = opt.largeNoiseSVPrior;
        logSmallNoisePrior = std::log(opt.smallNoiseSVPrior);
        logLargeNoisePrior = std::log(opt.largeNoiseSVPrior);
    }

    float
    logPrior(
        const unsigned gt,
        const float largeNoiseWeight) const
    {
        assert(largeNoiseWeight >= 0. && largeNoiseWeight <= 1.);

        if (gt != SOMATIC_GT::NOISE) return _logPrior[gt];

        if (largeNoiseWeight <= 0.) return logSmallNoisePrior;
        if (largeNoiseWeight >= 1.) return logLargeNoisePrior;

        return std::log((1-largeNoiseWeight)*smallNoisePrior + largeNoiseWeight*largeNoisePrior);
    }

private:
    std::array<float,SOMATIC_GT::SIZE> prior;
    std::array<float,SOMATIC_GT::SIZE> _logPrior;

    float smallNoisePrior;
    float largeNoisePrior;
    float logSmallNoisePrior;
    float logLargeNoisePrior;
};



struct SVScorer
{
    SVScorer(
        const GSCOptions& opt,
        const bam_header_info& header);

    /// gather supporting evidence and generate:
    /// 1) diploid quality score and genotype for SV candidate
    /// 2) somatic quality score
    void
    scoreSV(
        const SVCandidateSetData& svData,
        const std::vector<SVCandidateAssemblyData>& mjAssemblyData,
        const SVMultiJunctionCandidate& mjSV,
        const std::vector<bool>& isJunctionFiltered,
        const bool isSomatic,
        std::vector<SVModelScoreInfo>& mjModelScoreInfo,
        SVModelScoreInfo& mjJointModelScoreInfo,
        bool& isMJEvent);

    typedef std::shared_ptr<SVScorePairProcessor> pairProcPtr;
    typedef std::shared_ptr<bam_streamer> streamPtr;

private:

    void
    processExistingAltPairInfo(
        const PairOptions& pairOpt,
        const SVCandidateSetData& svData,
        const SVCandidate& sv,
        SVEvidence& evidence);

    void
    getSVAltPairSupport(
        const PairOptions& pairOpt,
        const SVCandidateAssemblyData& assemblyData,
        const SVCandidate& sv,
        SVEvidence& evidence,
        std::vector<pairProcPtr>& pairProcList);

    /// find spanning read support for the reference allele for sv candidate
    void
    getSVRefPairSupport(
        const PairOptions& pairOpt,
        const SVCandidate& sv,
        SVEvidence& evidence,
        std::vector<pairProcPtr>& pairProcList);

    /// find paired read support for ref and alt alleles
    void
    getSVPairSupport(
        const SVCandidateSetData& svData,
        const SVCandidateAssemblyData& assemblyData,
        const SVCandidate& sv,
        SVEvidence& evidence);

    /// find split read support for ref and alt alleles
    void
    getSVSplitReadSupport(
        const SVCandidateAssemblyData& assemblyData,
        const SVCandidate& sv,
        SVScoreInfo& ssInfo,
        SVEvidence& evidence);

    /// determine maximum depth and MQ0 frac in region around breakend of normal sample
    void
    getBreakendMaxMappedDepthAndMQ0(
        const bool isMaxDepth,
        const double cutoffDepth,
        const SVBreakend& bp,
        unsigned& maxDepth,
        float& MQ0Frac);

    /// apply all scoring models relevant to this event:
    void
    computeAllScoreModels(
        const bool isSomatic,
        const std::vector<JunctionCallInfo>& junctionData,
        SVModelScoreInfo& modelScoreInfo);

    /// shared information gathering steps of all scoring models
    void
    scoreSV(
        const SVCandidateSetData& svData,
        const SVCandidateAssemblyData& assemblyData,
        const SVCandidate& sv,
        SVScoreInfo& ssInfo,
        SVEvidence& evidence);


    const std::vector<bool> _isAlignmentTumor;
    const bool _isRNA;
    const CallOptionsShared _callOpt;
    const CallOptionsSharedDeriv _callDopt;
    const CallOptionsDiploid _diploidOpt;
    const CallOptionsDiploidDeriv _diploidDopt;
    const ReadScannerOptions _scanOpt;
    const SVRefinerOptions _refineOpt;
    const CallOptionsSomatic _somaticOpt;
    const CallOptionsSomaticDeriv _somaticDopt;
    const ChromDepthFilterUtil _dFilterDiploid;
    const ChromDepthFilterUtil _dFilterSomatic;
    SVLocusScanner _readScanner;

    std::vector<streamPtr> _bamStreams;
};
