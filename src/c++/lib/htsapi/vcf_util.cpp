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

#include "htsapi/vcf_util.hpp"

#include <cassert>
#include <cctype>
#include <ctime>
#include "blt_util/blt_exception.hpp"

#include <iostream>
#include <sstream>

std::ostream& vcf_fileDate(std::ostream& os)
{
  const time_t t(time(nullptr));
  struct tm*   ct(localtime(&t));
  assert(nullptr != ct);

  static const unsigned dsize(64);
  char                  datebuf[dsize];
  const size_t          ret(strftime(datebuf, dsize, "%Y%m%d", ct));
  assert(ret != 0);
  return os << datebuf;
}

void write_vcf_filter(std::ostream& os, const char* id, const char* desc)
{
  os << "##FILTER=<ID=" << id << ",Description=\"" << desc << "\">\n";
}

struct gt_parse_helper {
  // return is_valid_genotype
  static bool start(const char* gt, std::vector<int>& gti, const bool is_badend)
  {
    gti.clear();
    if (isdigit(*gt)) return digit(gt, gti, is_badend);

    switch (*gt) {
    case '.':
      return unknown(gt, gti, is_badend);
    default:
      return false;
    }
  }

private:
  static bool unknown(const char* gt, std::vector<int>& gti, const bool is_badend)
  {
    gt++;
    gti.push_back(-1);
    switch (*gt) {
    case '\0':
      return true;
    case '|':
    case '/':
      return sep(gt, gti, is_badend);
    default:
      return is_badend;
    }
  }

  static bool sep(const char* gt, std::vector<int>& gti, const bool is_badend)
  {
    gt++;
    if (isdigit(*gt)) return digit(gt, gti, is_badend);
    switch (*gt) {
    case '.':
      return unknown(gt, gti, is_badend);
    default:
      return false;
    }
  }

  static bool digit(const char* gt, std::vector<int>& gti, const bool is_badend)
  {
    int val(0);
    while (isdigit(*gt)) {
      val = val * 10 + static_cast<int>(*gt - '0');
      gt++;
    }
    gti.push_back(val);

    switch (*gt) {
    case '\0':
      return true;
    case '|':
    case '/':
      return sep(gt, gti, is_badend);
    default:
      return is_badend;
    }
  }
};

void parse_gt(const char* gt, std::vector<int>& gti, const bool is_allow_bad_end_char)
{
  if (!gt_parse_helper::start(gt, gti, is_allow_bad_end_char)) {
    std::ostringstream oss;
    oss << "Can't parse genotype string: '" << gt << "'";
    throw blt_exception(oss.str().c_str());
  }
}
