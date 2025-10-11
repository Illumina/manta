//
// Manta - Structural Variant and Indel Caller
// Copyright (c) 2013-2025 Illumina, Inc.
//
// This program is licensed under the terms of the Polyform strict license
//
// ***As far as the law allows, the software comes as is, without
// any warranty or condition, and the licensor will not be liable
// to you for any damages arising out of these terms or the use
// or nature of the software, under any kind of legal claim.***
//
// You should have received a copy of the PolyForm Strict License 1.0.0
// along with this program.  If not, see <https://polyformproject.org/licenses/strict/1.0.0>.
//
//

/// \file
/// \author Chris Saunders
///

#include "blt_util/stat_util.hpp"

#include "boost/math/distributions/chi_squared.hpp"

bool is_chi_sqr_reject(const double xsq, const unsigned df, const double alpha)
{
  assert(xsq >= 0);
  assert(df > 0);

  boost::math::chi_squared dist(df);
  return ((1. - boost::math::cdf(dist, xsq)) < alpha);

#if 0
    // alternate implementation (is one faster?):
    const double xsq_crit_val(boost::math::quantile(dist,1.-alpha));
    return xsq>xsq_crit_val;
#endif
}

bool is_lrt_reject_null(
    const double null_loghood, const double alt_loghood, const unsigned df, const double alpha)
{
  if (df == 0) return false;
  if (null_loghood > alt_loghood) return false;

  const double log_lrt(-2. * (null_loghood - alt_loghood));

  return is_chi_sqr_reject(log_lrt, df, alpha);
}
