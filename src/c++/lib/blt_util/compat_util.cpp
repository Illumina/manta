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

#include "compat_util.hpp"

#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <cstring>

#include <iostream>

double compat_round(const double x)
{
  if (x >= 0.) {
    return std::floor(x + 0.5);
  } else {
    return std::ceil(x - 0.5);
  }
}

const char* compat_basename(const char* str)
{
#ifdef _MSC_VER
  static const char pathsep('\\');
#else
  static const char pathsep('/');
#endif
  const char* res(strrchr(str, pathsep));
  if (nullptr == res) return str;
  return res + 1;
}
