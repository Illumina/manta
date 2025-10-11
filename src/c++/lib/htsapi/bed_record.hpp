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

#include "blt_util/seq_util.hpp"

#include <iosfwd>
#include <string>
#include <vector>

struct bed_record {
  bed_record() { clear(); }

  /// \brief Construct bed_record from bed file line contained in string s
  ///
  /// \return false on error
  ///
  bool set(const char* s);

  void clear()
  {
    chrom.clear();
    begin = 0;
    end   = 0;
    line  = nullptr;
  }

  /// Test if bed records are valid, where valid means having a strictly positive region size
  bool is_valid() const { return (begin < end); }

  std::string chrom;
  int         begin = 0;
  int         end   = 0;
  const char* line  = nullptr;
};

std::ostream& operator<<(std::ostream& os, const bed_record& bedr);
