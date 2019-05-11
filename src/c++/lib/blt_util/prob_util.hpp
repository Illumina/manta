//
// Manta - Structural Variant and Indel Caller
// Copyright (c) 2013-2019 Illumina, Inc.
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
/// \author Chris Saunders
///

#pragma once

#include <cassert>
#include <cmath>

#include <iterator>
#include <type_traits>

/// given a value on [-inf,inf], transform it
/// to a value on [0,1]
///
/// This uses the simple binomial version of the 'soft-max'
/// multinomial transformation:
///
/// out = exp(in) / (1 + exp(in))
///
/// TODO: probably some roundoff/precision changes that would improve handling of certain input ranges
///
template <typename FloatType>
FloatType softMaxTransform(const FloatType real)
{
  static_assert(std::is_floating_point<FloatType>::value, "Requires floating point type.");

  if (real <= 0.) {
    const FloatType er(std::exp(real));
    return (er / (1. + er));
  } else {
    const FloatType enr(std::exp(-real));
    return (1. / (enr + 1.));
  }
}

template <typename FloatType>
FloatType softMaxInverseTransform(const FloatType ranged)
{
  static_assert(std::is_floating_point<FloatType>::value, "Requires floating point type.");

  assert((ranged >= 0) && (ranged <= 1));
  if (ranged <= 0.5) {
    return std::log(ranged / (1 - ranged));
  } else {
    return -std::log((1 - ranged) / ranged);
  }
}

/// helper function for softMaxTransform to map [0,1] output to [rangedMin,rangedMax] range instead:
template <typename FloatType>
FloatType softMaxTransform(const FloatType real, const FloatType rangedMin, const FloatType rangedMax)
{
  const FloatType range(rangedMax - rangedMin);
  const FloatType offset(rangedMin);
  return softMaxTransform(real) * range + offset;
}

template <typename FloatType>
FloatType softMaxInverseTransform(
    const FloatType ranged, const FloatType rangedMin, const FloatType rangedMax)
{
  const FloatType range(rangedMax - rangedMin);
  const FloatType offset(rangedMin);
  return softMaxInverseTransform((ranged - offset) / range);
}

/// Given a set of N values on [-inf,inf], transform them to
/// the probability simplex sum(N+1) < 1, where the N+1 term
/// is implicit and does not represent a new degree of freedom
/// in the prob distro. The function returns the final value of
/// the probability distro.
///
/// this uses basic soft-max multinomial transformation:
///
/// out_i = exp(in_i) / (1 + sum{i}{exp(in_i)})
///
/// TODO: probably some roundoff/precision changes that would improve handling of certain input ranges
///
template <typename IterType>
typename std::iterator_traits<IterType>::value_type softMaxRangeTransform(
    const IterType begin, const IterType end)
{
  typedef typename std::iterator_traits<IterType>::value_type value_type;
  static_assert(std::is_floating_point<value_type>::value, "Requires iterator on floating point type.");

  value_type sum(1);
  for (IterType i(begin); i != end; ++i) {
    *i = std::exp(*i);
    sum += *i;
  }

  const value_type norm(1 / sum);
  value_type       sum2(0);
  for (IterType i(begin); i != end; ++i) {
    *i *= norm;
    sum2 += *i;
  }
  return (1. - sum2);
}

template <typename IterType>
void softMaxInverseRangeTransform(const IterType begin, const IterType end)
{
  typedef typename std::iterator_traits<IterType>::value_type value_type;
  static_assert(std::is_floating_point<value_type>::value, "Requires iterator on floating point type.");

  value_type sum(0);
  for (IterType i(begin); i != end; ++i) {
    assert((*i >= 0) && (*i <= 1));
    sum += *i;
  }
  assert((sum >= 0) && (sum <= 1));

  const value_type norm(1. / (1 - sum));
  for (IterType i(begin); i != end; ++i) {
    *i = std::log(*i * norm);
  }
}

/// Find more accurate complement of probability distro:
///
/// This function is setup assuming that 1 - prob[cgt] could
/// create significant loss of precision due to floating point
/// artifact, so we sum all prob[! cgt] instead. Typically this
/// becomes valuable as prob[cgt] approaches 1.
///
template <typename It>
typename std::iterator_traits<It>::value_type prob_comp(It begin, const It end, const unsigned cgt)
{
  typedef typename std::iterator_traits<It>::value_type FloatType;
  static_assert(std::is_floating_point<FloatType>::value, "Requires iterator on floating point type.");

  unsigned  i(0);
  FloatType val(0.);
  for (; begin != end; ++begin, ++i) {
    if (i == cgt) continue;
    val = val + *begin;
  }
  return val;
}

/// Normalize log-transformed probability distro
///
/// \param[out] maxElementIndex Index of the most probable component
///
template <typename It>
void normalizeLogDistro(const It pbegin, const It pend, unsigned& maxElementIndex)
{
  typedef typename std::iterator_traits<It>::value_type FloatType;
  static_assert(std::is_floating_point<FloatType>::value, "Requires iterator on floating point type.");

  // scale and exp pprob values:
  maxElementIndex = 0;
  if (pbegin == pend) return;
  FloatType max(*pbegin);
  unsigned  i(1);
  for (It p(pbegin + 1); p != pend; ++p, ++i) {
    if (*p > max) {
      max             = *p;
      maxElementIndex = i;
    }
  }

  FloatType sum(0.);
  for (It p(pbegin); p != pend; ++p) {
    *p = std::exp(*p - max);  // To alleviate underflow problem
    sum += *p;
  }

  // normalize:
  sum = 1. / sum;
  for (It p(pbegin); p != pend; ++p) {
    *p *= sum;
  }
}

/// optimized version of probability distro normalization for log-transformed input
///
/// optimization is to treat values significantly less than opt-max as zero probability
///
/// opt-max is found within the subset of the distribution where the predicate
/// iterator is true
///
/// \param[out] maxElementIndex Index of the most probable component
template <typename It, typename It2>
void opt_normalize_ln_distro(const It pbegin, const It pend, const It2 pred_begin, unsigned& maxElementIndex)
{
  typedef typename std::iterator_traits<It>::value_type FloatType;
  static_assert(std::is_floating_point<FloatType>::value, "Requires iterator on floating point type.");

  maxElementIndex = 0;
  if (pbegin == pend) return;

  bool      is_max(false), is_opt_max(false);
  FloatType max(0), opt_max(0);

  unsigned i(0);
  It2      pred(pred_begin);
  for (It p(pbegin); p != pend; ++p, ++pred, ++i) {
    if ((!is_max) || (*p > max)) {
      max             = *p;
      maxElementIndex = i;
      is_max          = true;
    }
    if (((!is_opt_max) || (*p > max)) && *pred) {
      opt_max    = *p;
      is_opt_max = true;
    }
  }

  assert(is_opt_max);

  static const FloatType norm_thresh(20);
  static const FloatType opt_thresh(5);

  FloatType sum(0.);
  pred = (pred_begin);
  for (It p(pbegin); p != pend; ++p, ++pred) {
    FloatType  mdiff(max - *p);
    const bool is_mdiff_skip(mdiff > norm_thresh);
    if (is_mdiff_skip) {
      if (!*pred) {
        *p = 0;
        continue;
      }
      FloatType optdiff(opt_max - *p);
      if (optdiff > opt_thresh) {
        *p = 0;
        continue;
      }
    }
    *p = std::exp(-mdiff);
    sum += *p;
  }

  // normalize:
  sum = 1. / sum;
  for (It p(pbegin); p != pend; ++p) {
    *p *= sum;
  }
}

void check_ln_distro_invalid_value(const char* label, const double val, const unsigned n);

void check_ln_distro_invalid_sum(const char* label, const double sum);

template <typename It>
double check_ln_distro(
    It i, const It i_end, const char* label, const double tol = 0.00001, const double target = 1)
{
  unsigned n(1);
  double   sum(0);
  for (; i != i_end; ++i, ++n) {
    const double val(std::exp(*i));
    if ((val < 0.) || (val > 1.)) {
      check_ln_distro_invalid_value(label, val, n);
    }
    sum += val;
  }
  if (std::abs(sum - target) > tol) {
    check_ln_distro_invalid_sum(label, sum);
  }
  return sum;
}
