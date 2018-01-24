//
// Manta - Structural Variant and Indel Caller
// Copyright (c) 2013-2017 Illumina, Inc.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//

/// \file
/// \author Chris Saunders and Xiaoyu Chen
///

#pragma once

#include <cassert>
#include <cmath>
#include <cstdlib>

#include <iosfwd>
#include <set>
#include <string>
#include <vector>


namespace DIPLOID_GT
{
enum index_t
{
    REF,
    HET,
    HOM,
    SIZE
};

static const float altPriors[] = {0., 0.5, 0.99};

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
        return nullptr;
    }
}

inline
const char*
label(const unsigned i)
{
    return label(static_cast<index_t>(i));
}

inline
float
altFraction(const index_t i)
{
    switch (i)
    {
    case REF :
        return altPriors[0];
    case HET :
        return altPriors[1];
    case HOM :
        return altPriors[2];
    default:
        assert(false && "Unknown GT state");
        return 0;
    }
}

inline
double
altLnFraction(const index_t i)
{
    switch (i)
    {
    case REF :
        return std::log(altPriors[0]);
    case HET :
        return std::log(altPriors[1]);
    case HOM :
        return std::log(altPriors[2]);
    default:
        assert(false && "Unknown GT state");
        return 0;
    }
}

inline
double
altLnCompFraction(const index_t i)
{
    switch (i)
    {
    case REF :
        return std::log(1-altPriors[0]);
    case HET :
        return std::log(1-altPriors[1]);
    case HOM :
        return std::log(1-altPriors[2]);
    default:
        assert(false && "Unknown GT state");
        return 0;
    }
}

}



struct SVScoreInfoDiploidSample
{
    SVScoreInfoDiploidSample() :
        phredLoghood(DIPLOID_GT::SIZE,0),
        pprob(DIPLOID_GT::SIZE,0)
    {}

    void
    clear()
    {
        filters.clear();
        gt=DIPLOID_GT::REF;
        gtScore=0;
        std::fill(phredLoghood.begin(),phredLoghood.end(),0);
        std::fill(pprob.begin(), pprob.end(),0);
    }

    std::set<std::string> filters;

    DIPLOID_GT::index_t gt = DIPLOID_GT::REF;

    unsigned gtScore = 0; ///< quality score of genotype

    std::vector<unsigned> phredLoghood;
    std::vector<double> pprob;
};


std::ostream&
operator<<(
    std::ostream& os,
    const SVScoreInfoDiploidSample& sid);



/// consolidate all germline scoring results applied to an SV candidate
struct SVScoreInfoDiploid
{
    void
    setSampleCount(
        const unsigned sampleCount)
    {
        samples.resize(sampleCount);
    }

    void
    clear()
    {
        filters.clear();
        altScore=0;
        for (auto& sample : samples)
        {
            sample.clear();
        }
    }

    std::set<std::string> filters;

    unsigned altScore = 0; ///< quality score indicating any non-reference state (regardless of specific genotype)

    std::vector<SVScoreInfoDiploidSample> samples;
};



std::ostream&
operator<<(
    std::ostream& os,
    const SVScoreInfoDiploid& sid);
