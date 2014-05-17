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
/// \author Chris Saunders and Xiaoyu Chen
///

#pragma once

#include "manta/SVCandidateAssemblyData.hh"
#include "manta/SVCandidate.hh"

#include <iosfwd>
#include <string>


struct SVAlignmentInfo
{
    SVAlignmentInfo(
        const SVCandidate& sv,
        const SVCandidateAssemblyData& assemblyData);

    bool
    isBp1SplitLeftOfHomologyRange() const
    {
        return _isBp1LeftOpen;
    }

    bool
    isBp2SplitLeftOfHomologyRange() const
    {
        return _isBp2LeftOpen;
    }

    bool
    isAltAlleleSpecificSequence() const
    {
        return _isAltAlleleSpecificSequence;
    }

    bool
    isRefAlleleSpecificSequence() const
    {
        return _isRefAlleleSpecificSequence;
    }

    const std::string&
    bp1ContigSeq() const
    {
        return (_isBp1ContigReversed ? revContigSeq : contigSeq);
    }

    const std::string&
    bp2ContigSeq() const
    {
        return (_isBp2ContigReversed ? revContigSeq : contigSeq);
    }

    const std::string&
    bp1ReferenceSeq() const
    {
        return bp1RefSeq;
    }

    const std::string&
    bp2ReferenceSeq() const
    {
        return (_isSpanning ? bp2RefSeq : bp1RefSeq);
    }

    bool
    isSpanning() const
    {
        return _isSpanning;
    }

    /// do we have a enough room on either side of the breakend for
    /// both the ref and contig to make a fair split read evaluation?
    bool
    isMinBpEdge(const unsigned minEdge) const;

    friend
    std::ostream&
    operator<<(std::ostream& os, const SVAlignmentInfo& ai);

private:
    std::string contigSeq;
    std::string revContigSeq;
    std::string bp1RefSeq;
    std::string bp2RefSeq;
    const bool _isSpanning;
    const bool _isBp1ContigReversed;
    const bool _isBp2ContigReversed;

    // breakpoints going through split read evaluation must be classified as left-open
    // or right-open only, so it's safe to capture only one bit for each breakend,
    // example of a deletion with bp1 occurring first in reference coordinate is
    // (right-open, left-open):
    const bool _isBp1LeftOpen;
    const bool _isBp2LeftOpen;

    // the presence of allele specific sequence (such as the deleted sequence in the
    // reference allele of a deletion) affects how we handle breakend homology for
    // that allele.
    //
    // note there are two bits here because in a complex indel both alleles can contain allele-specific sequence.
    //
    bool _isRefAlleleSpecificSequence;
    bool _isAltAlleleSpecificSequence;

public:
    /// all offset range 'begin' values correspond to the zero-indexed base immediately before the breakend on the fwd-strand,
    /// and 'end' values correspond to the zero-indexed base immediately before the breakend on the forward strand+microhomology range
    /// In the absence of microhomology, begin and end should be equal.
    known_pos_range2 bp1ContigOffset;
    known_pos_range2 bp2ContigOffset;
    known_pos_range2 bp1RefOffset;
    known_pos_range2 bp2RefOffset;
};

std::ostream&
operator<<(std::ostream& os, const SVAlignmentInfo& ai);



/// sample-specific and allele-specific evidence info
struct SVSampleAlleleInfo
{
    SVSampleAlleleInfo() :
        spanningPairCount(0),
        confidentSpanningPairCount(0),
        confidentSemiMappedSpanningPairCount(0),
        splitReadCount(0),
        splitReadEvidence(0),
        splitReadMapQ(0),
        confidentSplitReadCount(0)
    {}

    void
    clear()
    {
        spanningPairCount = 0;
        confidentSpanningPairCount = 0;
        confidentSemiMappedSpanningPairCount = 0;
        splitReadCount = 0;
        splitReadEvidence = 0;
        splitReadMapQ = 0;
        confidentSplitReadCount = 0;
    }

    // allele pair support
    unsigned spanningPairCount;  ///< all mapped pairs compatible with the allele
    unsigned confidentSpanningPairCount;  ///< pairs where both reads are mapped and we've successfully looked up a fragment prob of 0.01 or more
    unsigned confidentSemiMappedSpanningPairCount; ///< pairs where at least one read is mapped and we've successfully looked up a fragment prob of 0.01 or more

    // allele split support
    unsigned splitReadCount;
    float splitReadEvidence;
    float splitReadMapQ;

    unsigned confidentSplitReadCount; ///< count by comparing alignment quality vs the other allele
};

std::ostream&
operator<<(std::ostream& os, const SVSampleAlleleInfo& si);



/// sample-specific evidence info
struct SVSampleInfo
{
    void
    clear()
    {
        alt.clear();
        ref.clear();
    }

    SVSampleAlleleInfo alt;
    SVSampleAlleleInfo ref;
};

std::ostream&
operator<<(std::ostream& os, const SVSampleInfo& si);


/// consolidate model-agnostic scoring results applied to an SV candidate
struct SVScoreInfo
{
    SVScoreInfo() :
        bp1MaxDepth(0),
        bp2MaxDepth(0),
        bp1MQ0Frac(0),
        bp2MQ0Frac(0)
    {}

    void
    clear()
    {
        normal.clear();
        tumor.clear();

        bp1MaxDepth = 0;
        bp2MaxDepth = 0;

        bp1MQ0Frac = 0.;
        bp2MQ0Frac = 0.;
    }

    SVSampleInfo normal;
    SVSampleInfo tumor;

    unsigned bp1MaxDepth;
    unsigned bp2MaxDepth;

    float bp1MQ0Frac;
    float bp2MQ0Frac;
};


std::ostream&
operator<<(std::ostream& os, const SVScoreInfo& ssi);
