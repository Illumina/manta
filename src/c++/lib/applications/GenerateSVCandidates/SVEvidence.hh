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


///
/// \brief classes required to build the SVEvidence object
///
/// Note that the SVEvidence object is a kind of mini-database, accumulating all information about the
/// relationship of each fragment with a specific SV-candidate. This allows us to tie together:
///
/// (1) spanning read information from the full fragment
/// (2) split read information from read1
/// (3) split read information from read2
/// (4) fragment mapping properties from read1 & read2
///
/// ...for each fragment. The object stores this per-fragment information for all fragments impacting
/// all alleles at a given locus (But note at present we limit the alternate alleles to one).
///
/// The class structure below represents kind of a crummy db schema -- there's probably a better way to
/// do this?? Open to suggestions...
///




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
    float splitEvidence; ///< if evaluated, what is the evidence score?
};


/// track all support data from an individual fragment specific to an individual breakend of a single allele
///
struct SVFragmentEvidenceAlleleBreakend
{
    SVFragmentEvidenceAlleleBreakend() :
        isFragmentSupport(false),
        fragLengthProb(0)
    {}

    SVFragmentEvidenceAlleleBreakendPerRead&
    getRead(const bool isRead1)
    {
        return (isRead1 ? read1 : read2);
    }

    bool isFragmentSupport; ///< if true, paired-read analysis shows that this read pair fragment supports this allele on this breakend
    float fragLengthProb; ///< if isFragmentSupport, what is the prob of the fragment size given this allele?
    SVFragmentEvidenceAlleleBreakendPerRead read1; // read1 specific evidence
    SVFragmentEvidenceAlleleBreakendPerRead read2; // read2 specific evidence
};


/// track all support data from an individual fragment specific to a single allele of an SV candidate
///
struct SVFragmentEvidenceAllele
{
    SVFragmentEvidenceAlleleBreakend&
    getBp(const bool isBp1)
    {
        return (isBp1 ? bp1 : bp2 );
    }

    SVFragmentEvidenceAlleleBreakend bp1;
    SVFragmentEvidenceAlleleBreakend bp2;
};


/// store properties of the reads in a fragment which are not tightly coupled to any one allele/bp, etc....
///
struct SVFragmentEvidenceRead
{
    SVFragmentEvidenceRead() :
        isScanned(false),
        isAnchored(false),
        mapq(0)
    {}

    bool
    isObservedAnchor() const
    {
        return (isScanned && isAnchored);
    }

    bool isScanned; ///< if true, this read's bam record has been scanned to fill in the remaining values in this object

    bool isAnchored; ///< if true, the read is found and known to have a confident mapping wrt fragment support
    unsigned mapq;
};


/// track all support data from an individual fragment specific to an SV hypothesis
///
/// this is both to prevent double-counting of evidence and to consolidate different
/// sources of information (paired-split, etc).
///
struct SVFragmentEvidence
{
    SVFragmentEvidenceRead&
    getRead(const bool isRead1)
    {
        return (isRead1 ? read1 : read2);
    }

    SVFragmentEvidenceRead read1;
    SVFragmentEvidenceRead read2;

    SVFragmentEvidenceAllele alt;
    SVFragmentEvidenceAllele ref;
};


/// track all support data for an SV hypothesis
///
/// Note how this object is different than SVScoreInfoSomatic -- it is highly detailed and meant to be processed
/// to create summary statistics and scores later. Those scores and summary statistics should go into objects like
/// SomaticSVSCoreInfo to be written out in whichever output format is selected.
///
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
