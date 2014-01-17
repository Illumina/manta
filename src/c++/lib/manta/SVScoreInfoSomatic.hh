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
#include <cmath>


namespace SOMATIC_GT
{
	// TODO: estimated from tumor
	const double SOMATIC_MUTATION_FREQ = 0.6;

	enum index_t
	{
		REF,
		HET,
		HOM,
		SOM,
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
		case SOM :
			return "som";
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
	altFraction(const index_t i, const double somaticFreq)
	{
		switch (i)
		{
		case REF :
			return 0;
		case HET :
			return 0.5;
		case HOM :
			return 1.0;
		case SOM:
			//return SOMATIC_MUTATION_FREQ;
			return somaticFreq;

		default:
			assert(false && "Unknown GT state");
			return 0;
		}
	}

	inline
	double
	altLnFraction(const index_t i, const double somaticFreq)
	{
		//static const double val[] = { std::log(0.), std::log(0.5), std::log(1.), std::log(SOMATIC_MUTATION_FREQ) };
		static const double val[] = { std::log(0.), std::log(0.5), std::log(1.), std::log(somaticFreq) };
		switch (i)
		{
		case REF :
			return val[0];
		case HET :
			return val[1];
		case HOM :
			return val[2];
		case SOM:
			return val[3];
		default:
			assert(false && "Unknown GT state");
			return 0;
		}
	}

	inline
	double
	altLnCompFraction(const index_t i, const double somaticFreq)
	{
		return altLnFraction(static_cast<index_t>(2-i), somaticFreq);
	}
}


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

