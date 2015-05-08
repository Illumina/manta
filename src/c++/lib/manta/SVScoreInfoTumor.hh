// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Manta
// Copyright (c) 2013-2015 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/sequencing/licenses/>
//

///
/// \author Xiaoyu Chen
///

#pragma once

#include <cassert>
#include <cmath>
#include <cstdlib>

#include <iosfwd>
#include <set>
#include <string>


namespace TUMOR_GT
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
const char*
label(const unsigned i)
{
    return label(static_cast<const index_t>(i));
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
    	//TODO: fix this prior for low-AF
        return 0.5;
    case HOM :
        return 1.0;
    default:
        assert(false && "Unknown GT state");
        return 0;
    }
}

inline
double
altLnFraction(const index_t i)
{
	//TODO: fix value of HET for low-AF
    static const double val[] = { std::log(0.), std::log(0.5), std::log(1.) };
    switch (i)
    {
    case REF :
        return val[0];
    case HET :
        return val[1];
    case HOM :
        return val[2];
    default:
        assert(false && "Unknown GT state");
        return 0;
    }
}

inline
double
altLnCompFraction(const index_t i)
{
    return altLnFraction(static_cast<index_t>(2-i));
}

}

/// consolidate all tumor-only scoring results applied to an SV candidate
struct SVScoreInfoTumor
{
    void
    clear()
    {
        filters.clear();
        gt=TUMOR_GT::REF;
        altScore=0;
        gtScore=0;
    }

    std::set<std::string> filters;

    TUMOR_GT::index_t gt = TUMOR_GT::REF;

    unsigned altScore = 0; ///< quality score indicating any non-reference state (regardless of specific genotype)
    unsigned gtScore = 0; ///< quality score of genotype
};

std::ostream&
operator<<(
    std::ostream& os,
    const SVScoreInfoTumor& sid);
