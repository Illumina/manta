// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Manta
// Copyright (c) 2013 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/downloads/sequencing/licenses/>.
//

#include "boost/test/unit_test.hpp"

#include "parse_util.hh"

#include <string>

BOOST_AUTO_TEST_SUITE( parse_util )

using namespace illumina::blt_util;

BOOST_AUTO_TEST_CASE( test_parse_int ) {
    const char* two = "2";
    const int val(parse_int(two));
    BOOST_REQUIRE_EQUAL(val, 2);
}

BOOST_AUTO_TEST_CASE( test_parse_int_str ) {
    static const char two[] = "2";
    const int val(parse_int_str(std::string(two)));
    BOOST_REQUIRE_EQUAL(val, 2);
}

BOOST_AUTO_TEST_CASE( test_parse_int_str_bad_input ) {
    static const std::string junk("ABCD");
    BOOST_REQUIRE_THROW(parse_int_str(junk), std::exception);
}

BOOST_AUTO_TEST_SUITE_END()

