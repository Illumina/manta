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

#include "bed_record.hpp"
#include "blt_util/parse_util.hpp"

#include <cassert>

#include <iostream>

bool bed_record::set(const char* s)
{
  static const char     sep('\t');
  static const unsigned maxword(3);

  clear();

  line = s;

  // simple tab parse:
  const char* start(s);
  const char* p(start);

  unsigned wordindex(0);
  while (wordindex < maxword) {
    if ((*p == sep) || (*p == '\n') || (*p == '\0')) {
      switch (wordindex) {
      case 0:
        chrom = std::string(start, p - start);
        break;
      case 1:
        begin = illumina::blt_util::parse_int(start);
        assert(start == p);
        break;
      case 2:
        end = illumina::blt_util::parse_int(start);
        assert(start == p);
        break;
      default:
        assert(0);
        break;
      }
      start = p + 1;
      wordindex++;
    }
    if ((*p == '\n') || (*p == '\0')) break;
    ++p;
  }

  return (wordindex >= maxword);
}

std::ostream& operator<<(std::ostream& os, const bed_record& bedr)
{
  os << bedr.chrom << '\t' << bedr.begin << '\t' << bedr.end << '\n';

  return os;
}
