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

#include "blt_util/align_path.hpp"

//#define DEBUG_AP_TEST

#ifdef DEBUG_AP_TEST
#include <iostream>
namespace {
std::ostream& log_os(std::cerr);
}
#endif

BOOST_AUTO_TEST_SUITE(test_align_path)

using namespace ALIGNPATH;

static void test_single_cigar_conversion(const std::string& input)
{
  path_t apath;
  cigar_to_apath(input.c_str(), apath);
  BOOST_REQUIRE_EQUAL(input, apath_to_cigar(apath));
}

BOOST_AUTO_TEST_CASE(test_align_path_cigar_conversion)
{
  path_t apath;
  cigar_to_apath("10I10M10D10M10S", apath);
  BOOST_REQUIRE_EQUAL(apath.size(), 5u);

  // test round-trip:
  test_single_cigar_conversion("10I10M2S20M2I10M10D10M");
  test_single_cigar_conversion("");
  test_single_cigar_conversion("10S");
}

BOOST_AUTO_TEST_CASE(test_align_path_ref_length)
{
  path_t apath;
  cigar_to_apath("2I10M10D4I10M10N10M3S", apath);
  BOOST_REQUIRE_EQUAL(apath_ref_length(apath), 50u);
}

BOOST_AUTO_TEST_CASE(test_align_path_read_length)
{
  path_t apath;
  cigar_to_apath("2I10M10D4I10M10N10M3S", apath);
  BOOST_REQUIRE_EQUAL(apath_read_length(apath), 39u);
}

BOOST_AUTO_TEST_CASE(test_apath_indel_count)
{
  path_t apath;
  cigar_to_apath("2I10M10D4I10M1D10M10N10M3S", apath);
  BOOST_REQUIRE_EQUAL(apath_indel_count(apath), 3u);
}

static void test_string_clean(const char* cigar, const char* expect)
{
  path_t apath;
  cigar_to_apath(cigar, apath);
  apath_cleaner(apath);

  path_t expect_path;
  cigar_to_apath(expect, expect_path);

#ifdef DEBUG_AP_TEST
  log_os << "cleaned,expect: " << apath << " " << expect_path << "\n";
#endif

  BOOST_REQUIRE_EQUAL(apath, expect_path);
}

BOOST_AUTO_TEST_CASE(test_align_path_cleaner)
{
  // test cases for apath_cleaner function:
  test_string_clean("29M2I5M1D0M3I20M16S", "29M2I5M1D3I20M16S");
  test_string_clean("1M1P1M", "2M");
  test_string_clean("1M1D0M1I1M", "1M1D1I1M");
  //test_string_clean("0H1H1S0I1M1D0M1I1I1D1I1M1D1M","29M2I5M1D3I20M16S");
}

BOOST_AUTO_TEST_CASE(test_apath_clean_seqmatch)
{
  const std::string testCigar("10M1D10=2X10=1D1M1=1=1X1=1X");
  ALIGNPATH::path_t testPath;
  cigar_to_apath(testCigar.c_str(), testPath);

  apath_clean_seqmatch(testPath);

  BOOST_REQUIRE_EQUAL(apath_to_cigar(testPath), "10M1D22M1D6M");
}

BOOST_AUTO_TEST_CASE(test_apath_add_seqmatch)
{
  static const std::string testRead("AABAXXXY");
  static const std::string testRef("AAAADXXXX");

  static const std::string testCigar("4M1D4M");
  ALIGNPATH::path_t        testPath;
  cigar_to_apath(testCigar.c_str(), testPath);

  apath_add_seqmatch(testRead.begin(), testRead.end(), testRef.begin(), testRef.end(), testPath);

  BOOST_REQUIRE_EQUAL(apath_to_cigar(testPath), "2=1X1=1D3=1X");
}

static std::string test_limit_case(const std::string& cigar, const bool isReverse, const unsigned length)
{
  ALIGNPATH::path_t path;
  cigar_to_apath(cigar.c_str(), path);
  if (isReverse) {
    std::reverse(path.begin(), path.end());
  }
  apath_limit_ref_length(length, path);
  return apath_to_cigar(path);
}

BOOST_AUTO_TEST_CASE(test_apath_limit_ref_length)
{
  static const std::string testCigar("2=1X1=4D3=1X");
  BOOST_REQUIRE_EQUAL(test_limit_case(testCigar, false, 1), "1=");
  BOOST_REQUIRE_EQUAL(test_limit_case(testCigar, false, 5), "2=1X1=1D");
  BOOST_REQUIRE_EQUAL(test_limit_case(testCigar, false, 100), testCigar);

  BOOST_REQUIRE_EQUAL(test_limit_case(testCigar, true, 1), "1X");
  BOOST_REQUIRE_EQUAL(test_limit_case(testCigar, true, 5), "1X3=1D");

  static const std::string testCigar2("2=1X1=4I3=1X");
  BOOST_REQUIRE_EQUAL(test_limit_case(testCigar2, false, 5), "2=1X1=4I1=");
  BOOST_REQUIRE_EQUAL(test_limit_case(testCigar2, true, 5), "1X3=4I1=");
}

static std::string test_read_limit_case(
    const std::string& cigar, const bool isReverse, const unsigned start, const unsigned end)
{
  ALIGNPATH::path_t path;
  cigar_to_apath(cigar.c_str(), path);
  if (isReverse) {
    std::reverse(path.begin(), path.end());
  }
  apath_limit_read_length(start, end, path);
  return apath_to_cigar(path);
}

BOOST_AUTO_TEST_CASE(test_apath_limit_read_length)
{
  static const std::string testCigar("2=1X1=4I3=1X");
  BOOST_REQUIRE_EQUAL(test_read_limit_case(testCigar, false, 0, 1), "1=");
  BOOST_REQUIRE_EQUAL(test_read_limit_case(testCigar, false, 0, 5), "2=1X1=1I");
  BOOST_REQUIRE_EQUAL(test_read_limit_case(testCigar, false, 0, 100), testCigar);
  BOOST_REQUIRE_EQUAL(test_read_limit_case(testCigar, false, 1, 5), "1=1X1=1I");
  BOOST_REQUIRE_EQUAL(test_read_limit_case(testCigar, false, 2, 5), "1X1=1I");
  BOOST_REQUIRE_EQUAL(test_read_limit_case(testCigar, false, 5, 10), "3I2=");
}

BOOST_AUTO_TEST_CASE(test_apath_clip_trail)
{
  static const std::string testCigar("2S20M3S4H");
  ALIGNPATH::path_t        path;
  cigar_to_apath(testCigar.c_str(), path);

  BOOST_REQUIRE_EQUAL(apath_soft_clip_right_size(path), 3u);
}

BOOST_AUTO_TEST_SUITE_END()
