// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Manta
// Copyright (c) 2013-2014 Illumina, Inc.
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

#include "SVScorePairProcessor.hh"


struct SVScorePairAltProcessor : public SVScorePairProcessor
{
    SVScorePairAltProcessor(
        const std::vector<bool>& initIsAlignmentTumor,
        const SVLocusScanner& initReadScanner,
        const PairOptions& initPairOpt,
        const SVCandidate& initSv,
        const bool initIsBp1,
        SVEvidence& initEvidence) :
        SVScorePairProcessor(initIsAlignmentTumor, initReadScanner, initPairOpt, initSv, initIsBp1, initEvidence)
    {
        checkInput(sv);
    }

    void
    processClearedRecord(
        const bam_record& bamRead);

private:
    static
    void
    checkInput(
        const SVCandidate& sv);
};
