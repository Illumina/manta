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
