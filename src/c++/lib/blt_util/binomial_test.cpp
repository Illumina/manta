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

/// \file
/// \author Chris Saunders
///

#include "blt_util/binomial_test.hh"
#include "blt_util/stat_util.hh"

#include <boost/math/distributions/binomial.hpp>
#include <boost/math/distributions/complement.hpp>

using boost::math::binomial;
using boost::math::cdf;
using boost::math::complement;

#include <algorithm>



bool
is_reject_binomial_p_exact(const double alpha,
                           const double p,
                           const unsigned n_success,
                           const unsigned n_failure)
{

    const unsigned n_trial(n_success+n_failure);
    const double obs_p((double)n_success/(double)n_trial);

    double exact_prob;
    if (obs_p <= p)
    {
        exact_prob=cdf(binomial(n_trial,p),n_success);
    }
    else
    {
        exact_prob=cdf(binomial(n_trial,1.-p),n_failure);
    }

    return ((2.*exact_prob)<alpha);
}



bool
is_reject_binomial_p_chi_sqr(const double alpha,
                             const double p,
                             const unsigned n_success,
                             const unsigned n_failure)
{

    assert((p>0.) && (p<1.));

    const unsigned n_trial(n_success+n_failure);
    const double e_success(p*n_trial);
    const double e_failure(((double)n_trial)-e_success);

    const double d_success(n_success-e_success);
    const double d_failure(n_failure-e_failure);

    const double xsq((d_success*d_success)/e_success+(d_failure*d_failure)/e_failure);

    return is_chi_sqr_reject(xsq,1,alpha);
}

bool
is_reject_binomial_pval(const double threshold,
                        const double p,
                        const unsigned n_success,
                        const unsigned n_trials)
{
    assert(threshold >= 0);

    double observed_pval = get_binomial_pval(p, n_success, n_trials);

    return (observed_pval <= threshold);
}



double
get_binomial_pval(const double p,
                  const unsigned n_success,
                  const unsigned n_trials)
{
    //although technically binomial probabilities of
    // 0 or 1 are possible, they don't have much meaning
    assert((p >= 0.) && (p <= 1.));
    assert(n_success <= n_trials);
    if (n_success==0) return 1;

    return cdf(complement(binomial(n_trials, p), n_success - 1));
}



bool
is_reject_binomial_p(const double alpha,
                     const double p,
                     const unsigned n_success,
                     const unsigned n_failure)
{
    static const unsigned exact_test_threshold(250);

    const unsigned n_trial(n_success+n_failure);

    if (n_trial > exact_test_threshold)
    {
        return is_reject_binomial_p_chi_sqr(alpha,p,n_success,n_failure);
    }
    else
    {
        return is_reject_binomial_p_exact(alpha,p,n_success,n_failure);
    }
}
