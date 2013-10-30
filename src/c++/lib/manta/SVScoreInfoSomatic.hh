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

#include "manta/SVScoreInfo.hh"

#include <set>


/// consolidate all somatic scoring results applied to an SV candidate
struct SVScoreInfoSomatic : public SVScoreInfo
{
    typedef SVScoreInfo base_t;

    SVScoreInfoSomatic() :
        base_t(),
        somaticScore(0)
    {}

    void
    clear()
    {
        base_t::clear();

        filters.clear();
        somaticScore=0;
    }

    std::set<std::string> filters;

    unsigned somaticScore;
};


std::ostream&
operator<<(std::ostream& os, const SVScoreInfoSomatic& sis);

