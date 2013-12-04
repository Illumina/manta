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
/// \author Chris Saunders
///

#pragma once

#include "SVEvidence.hh"
#include "SVScorerPairOptions.hh"
#include "manta/BamRegionProcessor.hh"
#include "manta/SizeDistribution.hh"
#include "manta/SVCandidate.hh"
#include "manta/SVLocusScanner.hh"
#include "manta/SVScoreInfo.hh"

#include <vector>


struct SVScorePairAltInitParams
{
    SVScorePairAltInitParams(
        const SVLocusScanner& readScanner,
        const SVCandidate& sv,
        const bool isBp1);

    pos_t centerPos1;
    pos_t centerPos2;
    pos_t centerPos;

    // total impact of the alt allele on template size:
    pos_t altShift;
    unsigned minMapQ;
};


struct SVScorePairAltBamParams
{
    SVScorePairAltBamParams() :
        isSet(false),
        isTumor(false),
        minFrag(0),
        maxFrag(0),
        samplePtr(NULL),
        fragDistroPtr(NULL)
    {}

    bool isSet;
    bool isTumor;
    pos_t minFrag;
    pos_t maxFrag;
    SVSampleInfo* samplePtr;
    const SizeDistribution* fragDistroPtr;
    GenomeInterval interval;
};


struct SVScorePairAltProcessor : public BamRegionProcessor
{
    SVScorePairAltProcessor(
        const std::vector<bool>& initIsAlignmentTumor,
        const SVLocusScanner& initReadScanner,
        const PairOptions& initPairOpt,
        const SVCandidate& initSv,
        const bool initIsBp1,
        SVScoreInfo& initBaseInfo,
        SVEvidence& initEvidence);

    const GenomeInterval&
    nextBamIndex(
        const unsigned bamIndex);

    void
    processRecord(
        const bam_record& bamRead);

    const std::vector<bool> isAlignmentTumor;
    const SVLocusScanner& readScanner;
    const PairOptions& pairOpt;
    const SVCandidate& sv;
    const bool isBp1;
    SVScoreInfo& baseInfo;
    SVEvidence& evidence;

    const SVScorePairAltInitParams iparams;
    SVScorePairAltBamParams bparams;
};
