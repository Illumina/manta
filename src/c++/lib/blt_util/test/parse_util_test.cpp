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

#include "parse_util.hh"

#include <string>

BOOST_AUTO_TEST_SUITE( parse_util )

using namespace illumina::blt_util;


//
// check int parsing
//
BOOST_AUTO_TEST_CASE( test_parse_int )
{
    const char* two = "2";
    const int val(parse_int(two));
    BOOST_REQUIRE_EQUAL(val, 2);
}

BOOST_AUTO_TEST_CASE( test_parse_int_big )
{
    const char* twobig = "20000000000000000000";
    BOOST_REQUIRE_THROW(parse_int(twobig), std::exception);
}

BOOST_AUTO_TEST_CASE( test_parse_int_small )
{
    const char* twosmall = "-20000000000000000000";
    BOOST_REQUIRE_THROW(parse_int(twosmall), std::exception);
}

BOOST_AUTO_TEST_CASE( test_parse_int_empty )
{
    const char* empty = "";
    BOOST_REQUIRE_THROW(parse_int(empty), std::exception);
}

BOOST_AUTO_TEST_CASE( test_parse_int_tolerate_suffix )
{
    const char* suffix = "123abc";
    const int val(parse_int(suffix));
    BOOST_REQUIRE_EQUAL(val,123);
}

BOOST_AUTO_TEST_CASE( test_parse_int_str )
{
    static const char two[] = "2";
    const int val(parse_int_str(std::string(two)));
    BOOST_REQUIRE_EQUAL(val, 2);
}

BOOST_AUTO_TEST_CASE( test_parse_int_str_bad_input )
{
    static const std::string junk("ABCD");
    BOOST_REQUIRE_THROW(parse_int_str(junk), std::exception);
}


//
// check long parsing
//

BOOST_AUTO_TEST_CASE( test_parse_long )
{
    const char* two = "9223372036854775807";
    const long val(parse_long(two));
    BOOST_REQUIRE_EQUAL(val, 9223372036854775807l);
}

BOOST_AUTO_TEST_CASE( test_parse_long_big )
{
    const char* twobig = "9223372036854775808";
    BOOST_REQUIRE_THROW(parse_long(twobig), std::exception);
}

BOOST_AUTO_TEST_CASE( test_parse_long_small )
{
    const char* twosmall = "-20000000000000000000000000000000000000000000000000000";
    BOOST_REQUIRE_THROW(parse_long(twosmall), std::exception);
}

BOOST_AUTO_TEST_CASE( test_parse_long_empty )
{
    const char* empty = "";
    BOOST_REQUIRE_THROW(parse_long(empty), std::exception);
}

BOOST_AUTO_TEST_CASE( test_parse_long_str )
{
    static const char two[] = "2";
    const long val(parse_long_str(std::string(two)));
    BOOST_REQUIRE_EQUAL(val, 2l);
}


//
// check unsigned parsing
//
BOOST_AUTO_TEST_CASE( test_parse_unsigned )
{
    const char* two = "2";
    const unsigned val(parse_unsigned(two));
    BOOST_REQUIRE_EQUAL(val, 2u);
}

BOOST_AUTO_TEST_CASE( test_parse_unsigned_big )
{
    const char* twobig = "20000000000000000000";
    BOOST_REQUIRE_THROW(parse_unsigned(twobig), std::exception);
}

BOOST_AUTO_TEST_CASE( test_parse_unsigned_small )
{
    const char* twosmall = "-2";
    BOOST_REQUIRE_THROW(parse_unsigned(twosmall), std::exception);
}

BOOST_AUTO_TEST_CASE( test_parse_unsigned_empty )
{
    const char* empty = "";
    BOOST_REQUIRE_THROW(parse_unsigned(empty), std::exception);
}

BOOST_AUTO_TEST_CASE( test_parse_unsigned_tolerate_suffix )
{
    const char* suffix = "123abc";
    const unsigned val(parse_unsigned(suffix));
    BOOST_REQUIRE_EQUAL(val,123u);
}

BOOST_AUTO_TEST_CASE( test_parse_unsigned_str )
{
    static const char two[] = "2";
    const unsigned val(parse_unsigned_str(std::string(two)));
    BOOST_REQUIRE_EQUAL(val, 2u);
}

BOOST_AUTO_TEST_CASE( test_parse_unsigned_str_bad_input )
{
    static const std::string junk("ABCD");
    BOOST_REQUIRE_THROW(parse_unsigned_str(junk), std::exception);
}



//
// check double parsing
//
const double tol(0.0001);

BOOST_AUTO_TEST_CASE( test_parse_double )
{
    const char* two = "2.0";
    const double val(parse_double(two));
    BOOST_REQUIRE_CLOSE(val, 2.0, tol);
}

BOOST_AUTO_TEST_CASE( test_parse_double_exp )
{
    const char* big = "2.0e+100";
    const double val(parse_double(big));
    BOOST_REQUIRE_CLOSE(val, 2.0e+100, tol);
}

BOOST_AUTO_TEST_CASE( test_parse_double_inf )
{
    const char* inf = "Infinity";
    const double val(parse_double(inf));
    BOOST_REQUIRE(std::isinf(val));
}

BOOST_AUTO_TEST_CASE( test_parse_double_str )
{
    const char* two = "2.0";
    const double val(parse_double_str(std::string(two)));
    BOOST_REQUIRE_CLOSE(val, 2.0, tol);
}

BOOST_AUTO_TEST_CASE( test_parse_double_str_bad_input )
{
    static const std::string junk("ABCD");
    BOOST_REQUIRE_THROW(parse_double_str(junk), std::exception);
}

BOOST_AUTO_TEST_CASE( test_parse_double_str_bad_input2 )
{
    static const std::string junk("");
    BOOST_REQUIRE_THROW(parse_double_str(junk), std::exception);
}

BOOST_AUTO_TEST_SUITE_END()

