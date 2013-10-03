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

#include "blt_util/bam_streamer.hh"
#include "blt_util/bam_header_info.hh"
#include "manta/ChromDepthFilterUtil.hh"
#include "manta/SVCandidate.hh"
#include "manta/SVCandidateSetData.hh"
#include "manta/SVLocusScanner.hh"
#include "manta/SomaticSVScoreInfo.hh"
#include "assembly/AssembledContig.hh"
#include "manta/SVCandidateAssemblyData.hh"
#include "splitReadAlignment.hh"

#include "boost/shared_ptr.hpp"

#include <vector>
#include <string>


struct SVScorer
{
    SVScorer(
        const GSCOptions& opt,
        const bam_header_info& header);

    /// gather supporting evidence and generate somatic quality score for SV candidate
    void
    scoreSomaticSV(
        const SVCandidateSetData& svData,
        const SVCandidateAssemblyData& assemblyData,
        const SVCandidate& sv,
        SomaticSVScoreInfo& ssInfo);

private:

    /// find spanning read support for the reference allele in a single breakend
    void
    getSVRefPairSupport(
        const SVBreakend& bp,
        SomaticSVScoreInfo& ssInfo,
        const bool isBp1);

    /// find spanning read support for the reference allele for sv candidate
    void
    getSVRefPairSupport(
        const SVCandidate& sv,
        SomaticSVScoreInfo& ssInfo);

    /// find split read support for ref and alt alleles
    void
    getSVSplitReadSupport(
        const SVCandidateAssemblyData& assemblyData,
        const SVCandidate& sv,
        SomaticSVScoreInfo& ssInfo);

    /// determine maximum depth in region around breakend
    unsigned
    getBreakendMaxMappedDepth(
         const SVBreakend& bp);

    const std::vector<bool> _isAlignmentTumor;
    const SomaticCallOptions _somaticOpt;
    const ChromDepthFilterUtil _dFilter;
    SVLocusScanner _readScanner;

    typedef boost::shared_ptr<bam_streamer> streamPtr;
    std::vector<streamPtr> _bamStreams;
};

