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

#include "blt_util/compat_util.hpp"

#include <cassert>
#include <cstdlib>
#include <cstring>

#include <iosfwd>

struct ReadGroupLabel {
  /// if isCopyPtrs then the strings are copied and alloced/de-alloced by
  /// the object, if false the client is responsible these pointers over
  /// the lifetime of the label:
  ReadGroupLabel(const char* bamLabelInit, const char* rgLabelInit, const bool isCopyPtrsInit = true)
    : isCopyPtrs(isCopyPtrsInit),
      bamLabel((isCopyPtrs && (nullptr != bamLabelInit)) ? strdup(bamLabelInit) : bamLabelInit),
      rgLabel((isCopyPtrs && (nullptr != rgLabelInit)) ? strdup(rgLabelInit) : rgLabelInit)
  {
    assert(nullptr != bamLabel);
    assert(nullptr != rgLabel);
  }

  ReadGroupLabel(const ReadGroupLabel& rhs)
    : isCopyPtrs(rhs.isCopyPtrs),
      bamLabel(isCopyPtrs ? strdup(rhs.bamLabel) : rhs.bamLabel),
      rgLabel(isCopyPtrs ? strdup(rhs.rgLabel) : rhs.rgLabel)
  {
  }

  ReadGroupLabel& operator=(const ReadGroupLabel& rhs)
  {
    if (this == &rhs) return *this;
    clear();
    isCopyPtrs = rhs.isCopyPtrs;
    bamLabel   = (isCopyPtrs ? strdup(rhs.bamLabel) : rhs.bamLabel);
    rgLabel    = (isCopyPtrs ? strdup(rhs.rgLabel) : rhs.rgLabel);
    return *this;
  }

public:
  ~ReadGroupLabel() { clear(); }

  /// sort allowing for nullptr string pointers in primary and secondary key:
  bool operator<(const ReadGroupLabel& rhs) const
  {
    const int scval(strcmp(bamLabel, rhs.bamLabel));
    if (scval < 0) return true;
    if (scval == 0) {
      return (strcmp(rgLabel, rhs.rgLabel) < 0);
    }

    return false;
  }

private:
  void clear()
  {
    if (isCopyPtrs) {
      if (nullptr != bamLabel) free(const_cast<char*>(bamLabel));
      if (nullptr != rgLabel) free(const_cast<char*>(rgLabel));
    }
  }

  bool isCopyPtrs;

public:
  const char* bamLabel;
  const char* rgLabel;
};

std::ostream& operator<<(std::ostream& os, const ReadGroupLabel& rgl);
