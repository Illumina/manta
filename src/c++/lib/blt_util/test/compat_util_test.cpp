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

#include "boost/test/unit_test.hpp"

#include "compat_util.hpp"

#include <string>

BOOST_AUTO_TEST_SUITE(compat_util)

static void single_test_round(const double input, const double expect)
{
  static const double eps(0.00001);
  BOOST_CHECK_CLOSE(compat_round(input), expect, eps);
}

BOOST_AUTO_TEST_CASE(test_round)
{
  single_test_round(3.5, 4.0);
  single_test_round(3.2, 3.0);
  single_test_round(3.7, 4.0);
  single_test_round(-1.0, -1.0);
  single_test_round(-1.2, -1.0);
  single_test_round(-1.5, -2.0);
  single_test_round(-1.7, -2.0);
}

static void single_test_basename(const char* input, const char* expect)
{
  const char* result(compat_basename(input));
  BOOST_CHECK_EQUAL(std::string(result), std::string(expect));
}

BOOST_AUTO_TEST_CASE(test_basename)
{
  single_test_basename("foo", "foo");
  single_test_basename("", "");

#ifdef _WIN32
  single_test_basename("\\foo", "foo");
  single_test_basename("\\\\", "");
  single_test_basename("\\", "");
#else
  single_test_basename("/foo", "foo");
  single_test_basename("//", "");
  single_test_basename("/", "");
#endif
}

BOOST_AUTO_TEST_SUITE_END()
