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
