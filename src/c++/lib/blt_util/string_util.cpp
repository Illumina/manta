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

#include "string_util.hpp"

#include <cstring>

#include <iostream>

void split_string(const char* str, const char delimiter, std::vector<std::string>& v)
{
  v.clear();
  while (true) {
    const char* next(strchr(str, delimiter));
    if ((nullptr == next) || (delimiter == '\0')) {
      v.emplace_back(str);
      return;
    }
    v.emplace_back(str, next - str);
    str = next + 1;
  }
}

void destructive_split_string(char* str, const char delimiter, std::vector<const char*>& v)
{
  v.clear();
  while (true) {
    char* next(strchr(str, delimiter));
    v.push_back(str);
    if ((nullptr == next) || (delimiter == '\0')) return;
    *next = '\0';
    str   = next + 1;
  }
}

void split_string(
    const std::string& str, const char delimiter, std::vector<std::string>& v, const bool isSkipEmpty)
{
  v.clear();

  size_t start(0);
  while (true) {
    size_t next(str.find(delimiter, start));
    if (!(isSkipEmpty && ((next == start) || (next == std::string::npos)))) {
      v.emplace_back(str.substr(start, next - start));
    }
    if (next == std::string::npos) return;
    start = next + 1;
  }
}

bool split_match(const std::string& str, const char delimiter, const char* needle)
{
  size_t start(0);
  while (true) {
    size_t next(str.find(delimiter, start));
    if (0 == str.compare(start, next - start, needle)) return true;
    if (next == std::string::npos) break;
    start = next + 1;
  }
  return false;
}
