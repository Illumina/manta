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
#include "blt_util/ReadKey.hh"
#include "manta/ChromDepthFilterUtil.hh"
#include "manta/SVCandidate.hh"
#include "manta/SVCandidateSetData.hh"
#include "manta/SVLocusScanner.hh"
#include "manta/SVModelScoreInfo.hh"
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

    /// gather supporting evidence and generate:
    /// 1) diploid quality score and genotype for SV candidate
    /// 2) somatic quality score
    void
    scoreSV(
        const SVCandidateSetData& svData,
        const SVCandidateAssemblyData& assemblyData,
        const SVCandidate& sv,
        SVModelScoreInfo& modelScoreInfo);

private:

    /// track all support data from an individual fragment specific to a single allele of an SV hypothesis
    //
    struct SVFragmentEvidenceAllele
    {
        SVFragmentEvidenceAllele() :
            isFragmentSupport(false),
            FragmentSizeProb(0),
            isRead1SplitSupport(false),
            isRead2SplitSupport(false)
        {}

        bool isFragmentSupport; ///< if true, paired-read analysis shows that this read pair fragment supports this allele
        float FragmentSizeProb; ///< how likely is it to observe this fragment size, when fragment size is computed with this allele

        bool isRead1SplitSupport; ///< if true, read1 has been tested and shows support for this allele
        bool isRead2SplitSupport; ///< if true, read1 has been tested and shows support for this allele
    };

    /// track all support data from an individual fragment specific to an SV hypothesis
    //
    // this is both to prevent double-counting of evidnce and to consolidate different
    // sources of information (paired-split, etc).
    //
    struct SVFragmentEvidence
    {
        SVFragmentEvidenceAllele alt;
        SVFragmentEvidenceAllele ref;
    };

    // track all support data for an SV hypothesis
    //
    // Note how this object is different than SVScoreInfoSomatic -- it is highly detailed and meant to be processed
    // to create summary statistics and scores later. Those scores and summary statistics should go into objects like
    // SomaticSVSCoreInfo to be written out in whichever output format is selected.
    //
    struct SVEvidence
    {
        typedef std::map<std::string,SVFragmentEvidence> evidenceTrack_t;

        evidenceTrack_t normal;
        evidenceTrack_t tumor;
    };

    /// find spanning read support for the reference allele in a single breakend
    void
    getSVRefPairSupport(
        const SVBreakend& bp,
        SVScoreInfo& ssInfo,
        const bool isBp1);

    /// find spanning read support for the reference allele for sv candidate
    void
    getSVRefPairSupport(
        const SVCandidate& sv,
        SVScoreInfo& ssInfo);

    /// find paired read support for ref and alt alleles
    void
    getSVPairSupport(
        const SVCandidateSetData& svData,
        const SVCandidate& sv,
        SVScoreInfo& ssInfo);

    /// find split read support for ref and alt alleles
    void
    getSVSplitReadSupport(
        const SVCandidateAssemblyData& assemblyData,
        const SVCandidate& sv,
        SVScoreInfo& ssInfo);

    /// determine maximum depth in region around breakend
    unsigned
    getBreakendMaxMappedDepth(
        const SVBreakend& bp);


    /// shared information gathering steps of all scoring models
    void
    scoreSV(
        const SVCandidateSetData& svData,
        const SVCandidateAssemblyData& assemblyData,
        const SVCandidate& sv,
        SVScoreInfo& ssInfo,
        SVEvidence& evidence);


    const std::vector<bool> _isAlignmentTumor;
    const CallOptionsDiploid _diploidOpt;
    const CallOptionsSomatic _somaticOpt;
    const ChromDepthFilterUtil _dFilter;
    SVLocusScanner _readScanner;

    typedef boost::shared_ptr<bam_streamer> streamPtr;
    std::vector<streamPtr> _bamStreams;
};
