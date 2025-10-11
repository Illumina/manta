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

#include <map>

/// online median tracking obj assuming high repeat obs counts
///
/// Note that by design depth=0 is excluded from the median
struct MedianDepthTracker {
  void addObs(const unsigned val)
  {
    auto iter(_cmap.find(val));
    if (iter == _cmap.end()) {
      _cmap[val] = 1;
    } else {
      iter->second++;
    }
    _total++;
  }

  double getMedian() const
  {
    // +1 makes the 1/2 case work out correctly...
    unsigned   ztotal(_total + 1);
    const auto ziter(_cmap.find(0));
    if (ziter != _cmap.end()) {
      ztotal -= ziter->second;
    }

    unsigned sum        = 0;
    unsigned lastBefore = 0;
    unsigned firstAfter = 0;
    for (const auto& val : _cmap) {
      if (val.first == 0) continue;

      // double instead of half so that we stay away from float math:
      sum += (val.second * 2);
      if (sum >= ztotal) {
        firstAfter = val.first;
        if ((ztotal + val.second * 2) != (sum + 1)) {
          lastBefore = firstAfter;
        }
        break;
      }
      lastBefore = val.first;
    }

    assert((sum + 1) >= ztotal);

    if (lastBefore == firstAfter) return lastBefore;
    return (static_cast<double>(lastBefore + firstAfter) / 2.);
  }

private:
  unsigned                     _total = 0;
  std::map<unsigned, unsigned> _cmap;
};
