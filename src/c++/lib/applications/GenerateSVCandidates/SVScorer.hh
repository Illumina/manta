// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Manta
// Copyright (c) 2013 Illumina, Inc.
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
#include "SplitReadAlignment.hh"
#include "SVEvidence.hh"
#include "SVScorerPairOptions.hh"
#include "SVScorePairProcessor.hh"

#include "assembly/AssembledContig.hh"
#include "blt_util/bam_streamer.hh"
#include "blt_util/bam_header_info.hh"
#include "blt_util/qscore_snp.hh"
#include "manta/ChromDepthFilterUtil.hh"
#include "manta/SVCandidate.hh"
#include "manta/SVCandidateSetData.hh"
#include "manta/SVLocusScanner.hh"
#include "manta/SVModelScoreInfo.hh"
#include "manta/SVCandidateAssemblyData.hh"
#include "manta/SVScoreInfoSomatic.hh"

#include "boost/shared_ptr.hpp"

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


struct CallOptionsDiploidDeriv
{
    CallOptionsDiploidDeriv(
        const CallOptionsDiploid& opt)
    {
        using namespace DIPLOID_GT;

        assert(opt.indelPrior < 0.5);

        prior[HET] = opt.indelPrior;
        prior[HOM] = opt.indelPrior/2;
        prior[REF] = 1. - prior[HET] - prior[HOM];

        for (unsigned i(0); i<SIZE; ++i)
        {
            logPrior[i] = std::log(prior[i]);
        }
    }

    float prior[DIPLOID_GT::SIZE];
    float logPrior[DIPLOID_GT::SIZE];
};


struct CallOptionsSomaticDeriv
{
    CallOptionsSomaticDeriv(
        const CallOptionsSomatic& opt)
    {
        using namespace SOMATIC_GT;

        assert(opt.germlineSVPrior < 0.5);

        prior[SOM] = opt.somaticSVPrior;
        prior[HET] = opt.germlineSVPrior;
        prior[HOM] = opt.germlineSVPrior/2;
        prior[REF] = 1. - prior[SOM] - prior[HET] - prior[HOM];

        for (unsigned i(0); i<SIZE; ++i)
        {
            logPrior[i] = std::log(prior[i]);
        }
    }

    float prior[SOMATIC_GT::SIZE];
    float logPrior[SOMATIC_GT::SIZE];
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
        const SVCandidateAssemblyData& assemblyData,
        const SVCandidate& sv,
        const bool isSomatic,
        SVModelScoreInfo& modelScoreInfo);

    typedef boost::shared_ptr<SVScorePairProcessor> pairProcPtr;
    typedef boost::shared_ptr<bam_streamer> streamPtr;

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

    /// shared information gathering steps of all scoring models
    void
    scoreSV(
        const SVCandidateSetData& svData,
        const SVCandidateAssemblyData& assemblyData,
        const SVCandidate& sv,
        SVScoreInfo& ssInfo,
        SVEvidence& evidence);


    const std::vector<bool> _isAlignmentTumor;
    const CallOptionsShared _callOpt;
    const CallOptionsSharedDeriv _callDopt;
    const CallOptionsDiploid _diploidOpt;
    const CallOptionsDiploidDeriv _diploidDopt;
    const CallOptionsSomatic _somaticOpt;
    const CallOptionsSomaticDeriv _somaticDopt;
    const ChromDepthFilterUtil _dFilterDiploid;
    const ChromDepthFilterUtil _dFilterSomatic;
    SVLocusScanner _readScanner;

    std::vector<streamPtr> _bamStreams;
};
