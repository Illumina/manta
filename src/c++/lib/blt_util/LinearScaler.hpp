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

#include <cassert>

#include <algorithm>

/// \brief Maps an input to [0,1] based on where it lies within a [min,max] range specified during
/// initialization.
///
/// ### Example:
/// ```
/// LinearScalar<int> ls(50,100);
/// assert(ls.getScale(75) == 0.5);
/// ```
///
template <typename T>
struct LinearScaler {
  LinearScaler() : _min(static_cast<T>(0)), _factor(1.) {}

  LinearScaler(const T min, const T max) { init(min, max); }

  void init(const T min, const T max)
  {
    assert(max > min);
    _min    = min;
    _factor = (1. / static_cast<double>(max - min));
  }

  double getScale(const T val) const
  {
    static const double zero(0);
    static const double one(1);
    return std::min(one, std::max(zero, static_cast<double>(val - _min) * _factor));
  }

private:
  T      _min;
  double _factor;
};
