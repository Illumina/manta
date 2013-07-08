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
#include <vector>


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



/// consolidate all somatic scoring results applied to an SV candidate
struct SomaticSVScoreInfo
{
    SomaticSVScoreInfo() :
        somaticScore(0)
    {}

    void
    clear()
    {
        tumor.clear();
        normal.clear();
        filters.clear();
        somaticScore=0;
    }

    SVSampleInfo tumor;
    SVSampleInfo normal;
    std::vector<std::string> filters;

    unsigned somaticScore;
};


std::ostream&
operator<<(std::ostream& os, const SomaticSVScoreInfo& ssi);


