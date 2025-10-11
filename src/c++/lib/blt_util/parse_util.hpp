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

#include <string>

namespace illumina {
namespace blt_util {

/// parse c-string to TYPE
///
/// tolerates a non-TYPE suffix, but a non-empty prefix must be parsable as a TYPE,
/// on completion the value of s will reflect the extent of the parse
///
/// if available, specify s_end for minor performance improvement (in case of extremely large string)
///
unsigned parse_unsigned(const char*& s);

int parse_int(const char*& s);

long parse_long(const char*& s);

double parse_double(const char*& s, const char* s_end = nullptr);

/// parse c-string to TYPE
///
/// similar to above functions but:
/// - entire string must be convertible, no trailing suffix is allowed
/// - appropriate for rvalue char pointers
///
unsigned parse_unsigned_rvalue(const char* s);

int parse_int_rvalue(const char* s);

long parse_long_rvalue(const char* s);

double parse_double_rvalue(const char* s, const char* s_end = nullptr);

/// parse std::string to TYPE
///
/// entire string must be convertible, no trailing suffix is allowed
///
unsigned parse_unsigned_str(const std::string& s);

int parse_int_str(const std::string& s);

long parse_long_str(const std::string& s);

double parse_double_str(const std::string& s);

/// template version:
///
template <typename T>
T parse_type(const char*&)
{
  static_assert(sizeof(T) == 0, "no parse specialization available for type T");
  return T();
}

template <>
inline unsigned parse_type<unsigned>(const char*& s)
{
  return parse_unsigned(s);
}

template <>
inline int parse_type<int>(const char*& s)
{
  return parse_int(s);
}

template <>
inline long parse_type<long>(const char*& s)
{
  return parse_long(s);
}

template <>
inline double parse_type<double>(const char*& s)
{
  return parse_double(s);
}

}  // namespace blt_util
}  // namespace illumina
