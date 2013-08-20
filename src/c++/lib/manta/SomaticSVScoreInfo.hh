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
// <https://github.com/downloads/sequencing/licenses/>.
//

///
/// \author Chris Saunders
///

#pragma once

#include <iosfwd>
#include <string>
#include <set>


struct SVSampleInfo
{

    SVSampleInfo() :
        bp1SpanReads(0),
        bp2SpanReads(0),
        spanPairs(0)
    {}

    void
    clear()
    {
        bp1SpanReads=0;
        bp2SpanReads=0;
        spanPairs=0;
    }

    unsigned bp1SpanReads;
    unsigned bp2SpanReads;
    unsigned spanPairs;
};

std::ostream&
operator<<(std::ostream& os, const SVSampleInfo& si);

/// consolidate all somatic scoring results applied to an SV candidate
struct SomaticSVScoreInfo
{
    SomaticSVScoreInfo() :
        bp1MaxDepth(0),
        bp2MaxDepth(0),
        somaticScore(0)
    {}

    void
    clear()
    {
        tumor.clear();
        normal.clear();
        filters.clear();

        bp1MaxDepth=0;
        bp2MaxDepth=0;
        somaticScore=0;
    }

    SVSampleInfo tumor;
    SVSampleInfo normal;

    std::set<std::string> filters;

    unsigned bp1MaxDepth;
    unsigned bp2MaxDepth;
    unsigned somaticScore;
};


std::ostream&
operator<<(std::ostream& os, const SomaticSVScoreInfo& ssi);


