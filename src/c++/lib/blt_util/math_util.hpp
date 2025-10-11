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

#pragma once

#include "boost/math/special_functions/log1p.hpp"

#include <cmath>

#include <algorithm>

/// returns log(1+x), switches to special libc function when abs(x) is small
///
template <typename FloatType>
FloatType log1p_switch(const FloatType x)
{
  // better number??
  static const FloatType smallx_thresh(0.01);

  if (std::abs(x) < smallx_thresh) {
    return boost::math::log1p(x);
  } else {
    return std::log(1 + x);
  }
}

/// returns equiv of log(exp(x1)+exp(x2))
///
template <typename FloatType>
FloatType log_sum(FloatType x1, FloatType x2)
{
  if (x1 < x2) std::swap(x1, x2);
  return x1 + log1p_switch(std::exp(x2 - x1));
}

// helper for median() below
template <typename Iter>
typename std::iterator_traits<Iter>::value_type _ne_median(Iter begin, Iter end)
{
  assert(begin != end);
  const auto size(std::distance(begin, end));
  std::nth_element(begin, begin + size / 2, end);
  return *(begin + size / 2);
}

// helper for median() below
template <typename Iter>
typename std::iterator_traits<Iter>::value_type _ps_median(Iter begin, Iter end)
{
  assert(begin != end);
  const auto size(std::distance(begin, end));
  std::partial_sort(begin, begin + size / 2 + 1, end);
  return *(begin + size / 2);
}

/// returns median, partially reorders elements in specified range
///
template <typename Iter>
typename std::iterator_traits<Iter>::value_type median(Iter begin, Iter end)
{
  // Dispatch median call so as to avoid common broken libstdc++ impl

#ifndef BROKEN_NTH_ELEMENT
  // this is the preferred way to do it, it is optionally disabled because of common gcc bug:
  // https://gcc.gnu.org/bugzilla/show_bug.cgi?id=58800
  return _ne_median(begin, end);
#else
  return _ps_median(begin, end);
#endif
}

/// standardize the div by zero guard on division
///
/// note while this is a relatively simple utility, it was creeping
/// up all over the place with minor type variations. This standard
/// library copy should consistently do the sane thing for all integral
/// and floating point types, unless you want floating point wider than
/// double..
///
template <typename A, typename B>
double safeFrac(const A a, const B b)
{
  const double bd(static_cast<double>(b));
  return (((bd <= 0.) && (bd >= 0)) ? 0. : (a / bd));
}
