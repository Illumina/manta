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

#include "blt_util/blt_types.hpp"

#include <string>

/// \brief Manages a partial reference sequence segment
///
/// This object holds a subset of the reference sequence within a specific [begin,end] range,
/// plus some padding on each side. This scheme allows the client to store only the part of the
/// the reference that is currently required (to save memory), but access the reference using
/// the regular position coordinates of the full reference sequence.
///
/// \TODO Do not expose internal reference storage object type.
///
struct reference_contig_segment {
  reference_contig_segment() : _offset(0) {}

  char get_base(const pos_t pos) const
  {
    if (pos < _offset || pos >= end()) return 'N';
    return _seq[pos - _offset];
  }

  void get_substring(const pos_t pos, const pos_t length, std::string& substr) const
  {
    if (pos < _offset || (pos + length) > end()) {
      // slow path (minority of calls):
      substr.clear();
      for (int i(0); i < length; ++i) {
        substr.push_back(get_base(pos + i));
      }
    } else {
      // fast path
      substr.assign(_seq, pos - _offset, length);
    }
  }

  std::string&       seq() { return _seq; }
  const std::string& seq() const { return _seq; }

  pos_t get_offset() const { return _offset; }

  void set_offset(const pos_t offset) { _offset = offset; }

  pos_t end() const { return _offset + _seq.size(); }

  void clear()
  {
    _offset = 0;
    _seq.clear();
  }

private:
  pos_t       _offset;
  std::string _seq;
};
