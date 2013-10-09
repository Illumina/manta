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

#include <map>

/// track all support data from an individual read in a fragment specific to an individual breakend of a single allele
///
struct SVFragmentEvidenceAlleleBreakendPerRead
{
    SVFragmentEvidenceAlleleBreakendPerRead() :
        isSplitEvaluated(false),
        isSplitSupport(false),
        splitEvidence(0)
    {}

    bool isSplitEvaluated; ///< have we checked this read for split support of this bp?
    bool isSplitSupport; ///< if evaluated, does this read support this allele in the bp?
    float splitEvidence;
};

/// track all support data from an individual fragment specific to an individual breakend of a single allele
///
struct SVFragmentEvidenceAlleleBreakend
{
    SVFragmentEvidenceAlleleBreakend() :
        isFragmentSupport(false)
    {}

    SVFragmentEvidenceAlleleBreakendPerRead&
    getRead(const bool isRead1)
    {
        return (isRead1 ? read1 : read2);
    }

    bool isFragmentSupport; ///< if true, paired-read analysis shows that this read pair fragment supports this allele on this breakend
    SVFragmentEvidenceAlleleBreakendPerRead read1;
    SVFragmentEvidenceAlleleBreakendPerRead read2;
};


/// store properties intrinsic to the reads in a fragment
///
struct SVFragmentEvidenceRead
{
    SVFragmentEvidenceRead() :
        isScanned(false),
        isAnchored(false),
        mapq(0)
    {}

    bool isScanned; ///< if true, this read's bam record has been scanned to fill in the remaining values in this object

    bool isAnchored; ///< if true, the read is found and known to have a confident mapping wrt fragment support
    unsigned mapq;
};


/// track all support data from an individual fragment specific to a single allele of an SV hypothesis
//
struct SVFragmentEvidenceAllele
{
    SVFragmentEvidenceRead&
    getRead(const bool isRead1)
    {
        return (isRead1 ? read1 : read2);
    }

    SVFragmentEvidenceAlleleBreakend&
    getBp(const bool isBp1)
    {
        return (isBp1 ? bp1 : bp2 );
    }

    SVFragmentEvidenceRead read1;
    SVFragmentEvidenceRead read2;

    SVFragmentEvidenceAlleleBreakend bp1;
    SVFragmentEvidenceAlleleBreakend bp2;
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


    evidenceTrack_t&
    getSample(
        const bool isTumor)
    {
        return (isTumor ? tumor : normal);
    }

    evidenceTrack_t normal;
    evidenceTrack_t tumor;
};
