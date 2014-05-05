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

#include "manta/ShadowReadFinder.hh"
#include "manta/SVCandidateAssemblyData.hh"
#include "options/ReadScannerOptions.hh"
#include "options/SVRefinerOptions.hh"


struct ContigParams
{
    ContigParams(
        const SVCandidateAssemblyData& assemblyData,
        const SVCandidate& sv);

    // extended contig:
    const std::string& extSeq;

    /// where does the extended contig begin in reference coordinates:
    const int beginPos;

    /// where does the extended contig end in reference coordinates:
    const int endPos;

    known_pos_range2 bp1Offset;
    known_pos_range2 bp2Offset;
};


struct SVScorePairAltProcessor : public SVScorePairProcessor
{
    SVScorePairAltProcessor(
        const ReadScannerOptions& scanOpt,
        const SVRefinerOptions& refineOpt,
        const std::vector<bool>& initIsAlignmentTumor,
        const SVLocusScanner& initReadScanner,
        const PairOptions& initPairOpt,
        const SVCandidateAssemblyData& initAssemblyData,
        const SVCandidate& initSv,
        const bool initIsBp1,
        SVEvidence& initEvidence) :
        SVScorePairProcessor(initIsAlignmentTumor, initReadScanner, initPairOpt, initSv, initIsBp1, initEvidence),
        assemblyData(initAssemblyData),
        _shadowAligner(refineOpt.spanningAlignScores),
        _shadow(scanOpt.minSingletonMapqCandidates,
            (! initIsBp1), /// search for left-open shadows
            (  initIsBp1)), /// search for right-open shadows
        _contig(initAssemblyData,initSv)
    {
        checkInput(assemblyData,sv);
    }

    /// what to skip in addition to the core skip test?
    ///
    /// override to allow for shadow and chimera re-maps for large insertions:
    ///
    virtual
    bool
    isSkipRecord(
        const bam_record& bamRead)
    {
        if (! isLargeInsertSV()) return SVScorePairProcessor::isSkipRecord(bamRead);

        if (! bamRead.is_paired()) return true;
        else if (bamRead.is_unmapped() && bamRead.is_mate_unmapped()) return true;
        return false;
    }

    void
    processClearedRecord(
        const bam_record& bamRead);

private:
    static
    void
    checkInput(
        const SVCandidateAssemblyData& assemblyData,
        const SVCandidate& sv);

    static unsigned largeInsertSize() { return 100; }

    bool isLargeInsertSV() { return (sv.insertSeq.size() >= largeInsertSize() ); }


    bool
    alignShadowRead(
        const bam_record& bamRead,
        int& altTemplateSize);

    /// test whether a frag reference span provides sufficient support for a breakpoint of this sv:
    bool
    testFragOverlap(
        const int fragBeginRefPos,
        const int fragEndRefPos) const;

    ///////////////////////
    const SVCandidateAssemblyData& assemblyData;

    const GlobalAligner<int> _shadowAligner;
    ShadowReadFinder _shadow;

    ContigParams _contig;
};
