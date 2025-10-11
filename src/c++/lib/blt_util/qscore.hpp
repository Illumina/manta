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
/// \brief Various quality-score/probability conversion services
///
/// These methods rely on precomputed/cached conversion table to accelerate common quality score conversion
/// operations.
///
/// \author Chris Saunders
///

#pragma once

#include "blt_util/qscore_cache.hpp"

#include <cmath>

#include <algorithm>
#include <limits>
#include <type_traits>

template <typename FloatType>
FloatType error_prob_to_phred(const FloatType prob)
{
  static_assert(std::is_floating_point<FloatType>::value, "Transform requires floating point type.");
  static const FloatType minlog10(static_cast<FloatType>(std::numeric_limits<FloatType>::min_exponent10));
  return -10. * std::max(minlog10, std::log10(prob));
}

template <typename FloatType>
FloatType ln_error_prob_to_phred(const FloatType lnProb)
{
  static_assert(std::is_floating_point<FloatType>::value, "Transform requires floating point type.");
  static const FloatType minlog10(static_cast<FloatType>(std::numeric_limits<FloatType>::min_exponent10));
  static const FloatType ln10(std::log(static_cast<FloatType>(10)));
  return -10. * std::max(minlog10, lnProb / ln10);
}

template <typename FloatType>
int error_prob_to_qphred(const FloatType prob)
{
  return static_cast<int>(std::floor(error_prob_to_phred(prob) + 0.5));
}

template <typename FloatType>
int ln_error_prob_to_qphred(const FloatType lnProb)
{
  return static_cast<int>(std::floor(ln_error_prob_to_phred(lnProb) + 0.5));
}

inline double phred_to_error_prob(const double val)
{
  return std::pow(10., -val / 10.);
}

inline double qphred_to_error_prob(const int qscore)
{
  return qphred_cache::get_error_prob(qscore);
}

inline double qphred_to_ln_comp_error_prob(const int qscore)
{
  return qphred_cache::get_ln_comp_error_prob(qscore);
}

inline double qphred_to_ln_error_prob(const int qscore)
{
  return qphred_cache::get_ln_error_prob(qscore);
}

/// \brief Modify basecall error_prob score according to mapping quality of the read.
inline double phred_to_mapped_error_prob(const double basecall_val, const double mapping_val)
{
  const double be(phred_to_error_prob(basecall_val));
  const double me(phred_to_error_prob(mapping_val));
  return ((1. - me) * be) + (me * 0.75);
}

inline int qphred_to_mapped_qphred(const int basecall_val, const int mapping_val)
{
  return qphred_cache::get_mapped_qscore(basecall_val, mapping_val);
}
