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

/// \file

/// \author Chris Saunders
///

#pragma once

#include "blt_util/qscore_cache.hh"

#include <cmath>

#include <algorithm>
#include <limits>


template <typename FloatType>
int
error_prob_to_qphred(const FloatType prob)
{
    static const FloatType minlog10(static_cast<FloatType>(std::numeric_limits<FloatType>::min_exponent10));
    return static_cast<int>(std::floor(-10.*std::max(minlog10,std::log10(prob))+0.5));
}

inline
double
phred_to_error_prob(const double val)
{
    return std::pow(10.,-val/10.);
}

inline
double
qphred_to_error_prob(const int qscore)
{
    return qphred_cache::get_error_prob(qscore);
}

inline
double
qphred_to_ln_comp_error_prob(const int qscore)
{
    return qphred_cache::get_ln_comp_error_prob(qscore);
}

inline
double
qphred_to_ln_error_prob(const int qscore)
{
    return qphred_cache::get_ln_error_prob(qscore);
}


// modify basecall error_prob score according to mapping quality of the
// read:
inline
double
phred_to_mapped_error_prob(const double basecall_val,
                           const double mapping_val)
{
    const double be(phred_to_error_prob(basecall_val));
    const double me(phred_to_error_prob(mapping_val));
    return ((1.-me)*be)+(me*0.75);
}

inline
int
qphred_to_mapped_qphred(const int basecall_val,
                        const int mapping_val)
{
    return qphred_cache::get_mapped_qscore(basecall_val,mapping_val);
}


// logodds aka solexa:
inline
int
error_prob_to_qlogodds(const double prob)
{
    static const double maxlog10(static_cast<double>(std::numeric_limits<double>::max_exponent10));
    return static_cast<int>(std::floor(10.*std::min(maxlog10,std::log10((1.-prob)/prob))+0.5));
}

inline
double
logodds_to_error_prob(const double val)
{
    return 1./(1.+std::pow(10.,val/10.));
}

inline
double
qlogodds_to_error_prob(const int qscore)
{
    return qlogodds_cache::get_error_prob(qscore);
}

inline
int
error_prob_to_qscore(const double prob,
                     const bool is_qphred)
{
    if (is_qphred)
    {
        return error_prob_to_qphred(prob);
    }
    else
    {
        return error_prob_to_qlogodds(prob);
    }
}

inline
double
qscore_to_error_prob(const int qscore,
                     const bool is_qphred)
{
    if (is_qphred)
    {
        return qphred_to_error_prob(qscore);
    }
    else
    {
        return qlogodds_to_error_prob(qscore);
    }
}

