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

#include "blt_types.hpp"
#include "known_pos_range2.hpp"

#include "boost/icl/interval_map.hpp"

/// \brief Accumulates region specific sum(T) give a set of (region,T) associations.
///
/// Example:
/// Given the following set of input ranges and payload values:
/// range: [0,10] value: 3
/// range: [5,15] value: 1
/// range: [12,16] value: 2
///
/// This object internally constructs region-specific sum values:
/// range: [0,4] value: 3
/// range: [5,10] value: 4
/// range: [11,11] value: 1
/// range: [12,15] value: 3
/// range: [15,16] value: 2
///
/// ..and the maxVal() function below will return 4.
///
template <typename T>
struct RegionSum {
  void clear() { _map.clear(); }

  void add(const known_pos_range2 pr, const T val)
  {
    _map.add(std::make_pair(boost::icl::interval<pos_t>::right_open(pr.begin_pos(), pr.end_pos()), val));
  }

  /// return peak value for all regions:
  T maxVal() const
  {
    T    max(0);
    bool isFirst(true);
    for (const auto& val : _map) {
      if (isFirst || val.second > max) {
        max     = val.second;
        isFirst = false;
      }
    }
    return max;
  }

private:
  using map_t = boost::icl::interval_map<pos_t, T>;
  map_t _map;
};
