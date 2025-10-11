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

#include "blt_util/pos_range.hpp"

#include <iostream>

// output is always 1-indexed inclusive interval:
//
std::ostream& operator<<(std::ostream& os, const pos_range& pr)
{
  os << "[";
  if (pr.is_begin_pos) {
    os << pr.begin_pos + 1;
  } else {
    os << "-inf";
  }
  os << " .. ";
  if (pr.is_end_pos) {
    os << pr.end_pos;
  } else {
    os << "inf";
  }
  os << "]";

  return os;
}
