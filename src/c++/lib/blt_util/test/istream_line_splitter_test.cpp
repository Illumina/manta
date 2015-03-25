// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Manta
// Copyright (c) 2013-2015 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/sequencing/licenses/>
//

#include "boost/test/unit_test.hpp"

#include "istream_line_splitter.hh"

#include <sstream>
#include <string>


BOOST_AUTO_TEST_SUITE( test_istream_line_splitter )


BOOST_AUTO_TEST_CASE( test_istream_line_splitter_parse )
{

    std::string test_input("1\t2\t3\t4\n11\t22\t33\t44\n");
    std::istringstream iss(test_input);

    istream_line_splitter dparse(iss);

    int line_no(0);
    while (dparse.parse_line())
    {
        line_no++;
        static const unsigned expected_col_count(4);
        BOOST_CHECK_EQUAL(dparse.n_word(),expected_col_count);
        if       (1==line_no)
        {
            BOOST_CHECK_EQUAL(std::string(dparse.word[1]),std::string("2"));
        }
        else if (2==line_no)
        {
            BOOST_CHECK_EQUAL(std::string(dparse.word[1]),std::string("22"));
        }
    }
}


static
void
check_long_line(const int init_buffer_size)
{

    std::string test_input("1ABCDEFGHIJKLMNOPQRSTUVWXYZ\t2\t3\t4ABCDEFG\n11\t22\t33\t44XYZ\n");
    std::istringstream iss(test_input);

    istream_line_splitter dparse(iss,init_buffer_size);

    int line_no(0);
    while (dparse.parse_line())
    {
        line_no++;
        static const unsigned expected_col_count(4);
        BOOST_CHECK_EQUAL(dparse.n_word(),expected_col_count);
        if       (1==line_no)
        {
            BOOST_CHECK_EQUAL(std::string(dparse.word[0]),std::string("1ABCDEFGHIJKLMNOPQRSTUVWXYZ"));
            BOOST_CHECK_EQUAL(std::string(dparse.word[3]),std::string("4ABCDEFG"));
        }
        else if (2==line_no)
        {
            BOOST_CHECK_EQUAL(std::string(dparse.word[2]),std::string("33"));
        }
    }
}


BOOST_AUTO_TEST_CASE( test_istream_line_splitter_long_line )
{
    check_long_line(2);
    check_long_line(41);
}

BOOST_AUTO_TEST_SUITE_END()

