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


namespace DIPLOID_GT
{
enum index_t
{
    REF,
    HET,
    HOM,
    SIZE
};

inline
const char*
label(const index_t i)
{
    switch (i)
    {
    case REF :
        return "ref";
    case HET :
        return "het";
    case HOM :
        return "hom";
    default:
        assert(false && "Unknown GT state");
        return NULL;
    }
}

inline
float
altFraction(const index_t i)
{
    switch (i)
    {
    case REF :
        return 0;
    case HET :
        return 0.5;
    case HOM :
        return 1.0;
    default:
        assert(false && "Unknown GT state");
        return 0;
    }

}

inline
float
altFraction(const unsigned i)
{
    return altFraction(static_cast<const index_t>(i));
}

}


/// consolidate all germline scoring results applied to an SV candidate
struct SVScoreInfoDiploid : public SVScoreInfo
{
    typedef SVScoreInfo base_t;

    SVScoreInfoDiploid() :
        base_t(),
        gt(DIPLOID_GT::REF),
        altScore(0),
        gtScore(0)
    {}

    void
    clear()
    {
        base_t::clear();
        gtScore=0;
    }

    DIPLOID_GT::index_t gt;

    unsigned altScore; ///< quality score indicating any non-reference state (regardless of specific genotype)
    unsigned gtScore; ///< quality score of genotype
};

std::ostream&
operator<<(std::ostream& os, const SVScoreInfoDiploid& sid);
