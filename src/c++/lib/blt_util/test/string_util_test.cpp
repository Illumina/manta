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

#include "string_util.hh"

BOOST_AUTO_TEST_SUITE( string_util )

static const char* test_string("234562342");

template <typename T>
void
test_split_string_bytype(T test) {
    std::vector<std::string> result;

    split_string(test,'2',result);
    BOOST_CHECK_EQUAL(static_cast<int>(result.size()), 4);
    BOOST_CHECK_EQUAL(result[0], "");
    BOOST_CHECK_EQUAL(result[1], "3456");
    BOOST_CHECK_EQUAL(result[3], "");

    split_string(test,'X',result);
    BOOST_CHECK_EQUAL(static_cast<int>(result.size()), 1);
    BOOST_CHECK_EQUAL(result[0], test);

    split_string("",'X',result);
    BOOST_CHECK_EQUAL(static_cast<int>(result.size()), 1);
    BOOST_CHECK_EQUAL(result[0], "");
}

BOOST_AUTO_TEST_CASE( test_split_string_cstr ) {
    test_split_string_bytype(test_string);
}

BOOST_AUTO_TEST_CASE( test_split_string ) {
    const std::string test(test_string);
    test_split_string_bytype(test);
}

BOOST_AUTO_TEST_CASE( test_split_match ) {
    const std::string test(test_string);

    BOOST_CHECK(split_match(test,'2',"34"));
    BOOST_CHECK(! split_match(test,'2',"XX"));
}

BOOST_AUTO_TEST_SUITE_END()

