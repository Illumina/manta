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
