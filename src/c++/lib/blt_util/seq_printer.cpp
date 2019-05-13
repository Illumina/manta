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
