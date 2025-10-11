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

#include "boost/test/unit_test.hpp"

#include "istream_line_splitter.hpp"

#include <sstream>
#include <string>

BOOST_AUTO_TEST_SUITE(test_istream_line_splitter)

BOOST_AUTO_TEST_CASE(test_istream_line_splitter_parse)
{
  std::string        test_input("1\t2\t3\t4\n11\t22\t33\t44\n");
  std::istringstream iss(test_input);

  istream_line_splitter dparse(iss);

  int line_no(0);
  while (dparse.parse_line()) {
    line_no++;
    static const unsigned expected_col_count(4);
    BOOST_CHECK_EQUAL(dparse.n_word(), expected_col_count);
    if (1 == line_no) {
      BOOST_CHECK_EQUAL(std::string(dparse.word[1]), std::string("2"));
    } else if (2 == line_no) {
      BOOST_CHECK_EQUAL(std::string(dparse.word[1]), std::string("22"));
    }
  }
}

static void check_long_line(const int init_buffer_size)
{
  std::string        test_input("1ABCDEFGHIJKLMNOPQRSTUVWXYZ\t2\t3\t4ABCDEFG\n11\t22\t33\t44XYZ\n");
  std::istringstream iss(test_input);

  istream_line_splitter dparse(iss, init_buffer_size);

  int line_no(0);
  while (dparse.parse_line()) {
    line_no++;
    static const unsigned expected_col_count(4);
    BOOST_CHECK_EQUAL(dparse.n_word(), expected_col_count);
    if (1 == line_no) {
      BOOST_CHECK_EQUAL(std::string(dparse.word[0]), std::string("1ABCDEFGHIJKLMNOPQRSTUVWXYZ"));
      BOOST_CHECK_EQUAL(std::string(dparse.word[3]), std::string("4ABCDEFG"));
    } else if (2 == line_no) {
      BOOST_CHECK_EQUAL(std::string(dparse.word[2]), std::string("33"));
    }
  }
}

BOOST_AUTO_TEST_CASE(test_istream_line_splitter_long_line)
{
  check_long_line(2);
  check_long_line(41);
}

BOOST_AUTO_TEST_SUITE_END()
