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

#include "manta/SVCandidateAssemblyData.hh"
#include "manta/SVCandidate.hh"
#include "blt_util/ReadKey.hh"

#include <iosfwd>
#include <string>
#include <set>


struct SVAlignmentInfo
{
    SVAlignmentInfo(
        const SVCandidate& sv,
        const SVCandidateAssemblyData& assemblyData);

    const std::string&
    bp1ContigSeq() const
    {
        return (bp1ContigReversed ? revContigSeq : contigSeq);
    }

    const std::string&
    bp2ContigSeq() const
    {
        return (bp2ContigReversed ? revContigSeq : contigSeq);
    }

    friend
    std::ostream&
    operator<<(std::ostream& os, const SVAlignmentInfo& ai);

private:
    std::string contigSeq;
    std::string revContigSeq;
    bool bp1ContigReversed;
    bool bp2ContigReversed;

public:
    std::string bp1RefSeq;
    std::string bp2RefSeq;
    unsigned bp1ContigOffset;
    unsigned bp2ContigOffset;
    pos_t bp1RefOffset;
    pos_t bp2RefOffset;
};

std::ostream&
operator<<(std::ostream& os, const SVAlignmentInfo& ai);



/// sample-specific and allele-specific evidence info
struct SVSampleAlleleInfo
{
    SVSampleAlleleInfo() :
        bp1SpanReadCount(0),
        bp2SpanReadCount(0),
        spanPairCount(0),
        splitReadCount(0),
        splitReadEvidence(0),
        splitReadMapQ(0)
    {}

    void
    clear()
    {
        bp1SpanReadCount = 0;
        bp2SpanReadCount = 0;
        spanPairCount = 0;
        splitReadCount = 0;
        splitReadEvidence = 0;
        splitReadMapQ = 0;
    }

    // allele pair support
    unsigned bp1SpanReadCount;
    unsigned bp2SpanReadCount;
    unsigned spanPairCount;

    // allele split support
    unsigned splitReadCount;
    float splitReadEvidence;
    float splitReadMapQ;
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
        bp2MaxDepth(0)
    {}

    void
    clear()
    {
        normal.clear();
        tumor.clear();
        filters.clear();

        bp1MaxDepth=0;
        bp2MaxDepth=0;
    }

    SVSampleInfo normal;
    SVSampleInfo tumor;

    std::set<std::string> filters;

    unsigned bp1MaxDepth;
    unsigned bp2MaxDepth;
};


std::ostream&
operator<<(std::ostream& os, const SVScoreInfo& ssi);
