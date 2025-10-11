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

#include "blt_util/seq_printer.hpp"

#include <cassert>
#include <cstring>
#include <iostream>

void printSeq(const char* seq, std::ostream& os)
{
  static const unsigned rowSize(100);
  static const unsigned sectionSize(10);

  assert(nullptr != seq);
  const unsigned seqLen(strlen(seq));

  for (unsigned i(0); i < seqLen; ++i) {
    if (i) {
      if (0 == (i % rowSize))
        os << '\n';
      else if (0 == (i % sectionSize))
        os << ' ';
    }
    os << seq[i];
  }
}
