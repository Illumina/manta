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


struct SVScorePairInitParams
{
    SVScorePairInitParams(
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


struct SVScorePairBamParams
{
    SVScorePairBamParams() :
        isSet(false),
        isTumor(false),
        minFrag(0),
        maxFrag(0),
        fragDistroPtr(NULL)
    {}

    bool isSet;
    bool isTumor;
    pos_t minFrag;
    pos_t maxFrag;
    const SizeDistribution* fragDistroPtr;
    GenomeInterval interval;
};


struct SVScorePairProcessor : public BamRegionProcessor
{
    SVScorePairProcessor(
        const std::vector<bool>& initIsAlignmentTumor,
        const SVLocusScanner& initReadScanner,
        const PairOptions& initPairOpt,
        const SVCandidate& initSv,
        const bool initIsBp1,
        SVEvidence& initEvidence) :
        isAlignmentTumor(initIsAlignmentTumor),
        readScanner(initReadScanner),
        pairOpt(initPairOpt),
        sv(initSv),
        isBp1(initIsBp1),
        evidence(initEvidence),
        iparams(readScanner, sv, isBp1)
    {}

    const GenomeInterval&
    nextBamIndex(
        const unsigned bamIndex);

    void
    processRecord(
        const bam_record& bamRead)
    {
        if (isSkipRecord(bamRead)) return;
        processClearedRecord(bamRead);
    }

    // alternate interface
    static
    bool
    isSkipRecord(
        const bam_record& bamRead)
    {
        if (SVLocusScanner::isReadFilteredCore(bamRead)) return true;
        else if (bamRead.is_unmapped() || bamRead.is_mate_unmapped()) return true;
        else if (! is_innie_pair(bamRead)) return true;
        return false;
    }

    // process a record for which isSkipRecord() == false
    virtual
    void
    processClearedRecord(
        const bam_record& bamRead) = 0;

protected:

    static
    void
    setAlleleFrag(
        const SizeDistribution& fragDistro,
        const int size,
        SVFragmentEvidenceAlleleBreakend& bp)
    {
        float fragProb(fragDistro.cdf(size));
        fragProb = std::min(fragProb, (1-fragProb));
#ifdef DEBUG_MEGAPAIR
        log_os << __FUNCTION__ << ": fraglen,prob " << size << " " << fragProb << "\n";
#endif

        bp.isFragmentSupport = true;
        bp.fragLengthProb = fragProb;
    }

    const std::vector<bool> isAlignmentTumor;
    const SVLocusScanner& readScanner;
    const PairOptions& pairOpt;
    const SVCandidate& sv;
    const bool isBp1;
    SVEvidence& evidence;

    const SVScorePairInitParams iparams;
    SVScorePairBamParams bparams;
};
