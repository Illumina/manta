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
