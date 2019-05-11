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

#include "blt_util/parse_util.hpp"
#include "blt_util/blt_exception.hpp"

#include "boost/spirit/include/qi.hpp"

#include <cerrno>
#include <climits>
#include <cstdlib>

#include <limits>
#include <sstream>

static void parse_exception(const char* type_label, const char* parse_str)
{
  std::ostringstream oss;
  oss << "Can't parse " << type_label << " from string: '" << parse_str << "'";
  throw blt_exception(oss.str().c_str());
}

namespace illumina {
namespace blt_util {

static unsigned parse_unsigned_core(const char* s, const char*& s_out)
{
  static const int base(10);

  errno = 0;

  char*               endptr;
  const unsigned long val(strtoul(s, &endptr, base));
  if ((errno == ERANGE && (val == ULONG_MAX || val == 0)) || (errno != 0 && val == 0) || (endptr == s)) {
    parse_exception("unsigned long", s);
  }

  if (val > std::numeric_limits<unsigned>::max()) {
    parse_exception("unsigned", s);
  }

  s_out = endptr;

  return static_cast<unsigned>(val);
}

unsigned parse_unsigned(const char*& s)
{
  return parse_unsigned_core(s, s);
}

unsigned parse_unsigned_rvalue(const char* s)
{
  const char*    s_tmp(0);
  const unsigned val(parse_unsigned_core(s, s_tmp));
  if (*s_tmp != '\0') {
    parse_exception("unsigned", s);
  }
  return val;
}

unsigned parse_unsigned_str(const std::string& s)
{
  return parse_unsigned_rvalue(s.c_str());
}

static int parse_int_core(const char* s, const char*& s_out)
{
  const char* endptr(s);
  const long  val(parse_long(endptr));

  if ((val > std::numeric_limits<int>::max()) || (val < std::numeric_limits<int>::min())) {
    parse_exception("int", s);
  }

  s_out = endptr;

  return static_cast<int>(val);
}

int parse_int(const char*& s)
{
  return parse_int_core(s, s);
}

int parse_int_rvalue(const char* s)
{
  const char* s_tmp(0);
  const int   val(parse_int_core(s, s_tmp));
  if (*s_tmp != '\0') {
    parse_exception("int", s);
  }
  return val;
}

int parse_int_str(const std::string& s)
{
  return parse_int_rvalue(s.c_str());
}

static long parse_long_core(const char* s, const char*& s_out)
{
  static const int base(10);

  errno = 0;

  char*      endptr;
  const long val(strtol(s, &endptr, base));
  if ((errno == ERANGE && (val == LONG_MIN || val == LONG_MAX)) || (errno != 0 && val == 0) ||
      (endptr == s)) {
    parse_exception("long int", s);
  }

  s_out = endptr;

  return val;
}

long parse_long(const char*& s)
{
  return parse_long_core(s, s);
}

long parse_long_rvalue(const char* s)
{
  const char* s_tmp(0);
  const long  val(parse_long_core(s, s_tmp));
  if (*s_tmp != '\0') {
    parse_exception("long int", s);
  }
  return val;
}

long parse_long_str(const std::string& s)
{
  return parse_long_rvalue(s.c_str());
}

static double parse_double_core(const char* s, const char*& s_out, const char* s_end)
{
  double val;
  s_out = s;
  if (s_end == nullptr) s_end = s + strlen(s);
  bool isPass(boost::spirit::qi::parse(s_out, s_end, boost::spirit::double_, val));
  if (isPass) {
    isPass = (s != s_out);
  }
  if (!isPass) {
    parse_exception("double", s);
  }
  return val;
}

double parse_double(const char*& s, const char* s_end)
{
  return parse_double_core(s, s, s_end);
}

double parse_double_rvalue(const char* s, const char* s_end)
{
  const char*  s_tmp(0);
  const double val(parse_double_core(s, s_tmp, s_end));
  if (*s_tmp != '\0') {
    parse_exception("double", s);
  }
  return val;
}

double parse_double_str(const std::string& s)
{
  const char*       s2(s.c_str());
  const char* const s2_end(s2 + s.size());
  return parse_double_rvalue(s2, s2_end);
}

}  // namespace blt_util
}  // namespace illumina
