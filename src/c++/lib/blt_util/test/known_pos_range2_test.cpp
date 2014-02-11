// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Manta
// Copyright (c) 2013-2014 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/sequencing/licenses/>
//


#include "boost/test/unit_test.hpp"

#include "blt_util/known_pos_range2.hh"


BOOST_AUTO_TEST_SUITE( test_known_pos_range2 )


BOOST_AUTO_TEST_CASE( test_known_pos_range2_is_pos_intersect )
{
    // this corresponds to zero-index range [9,19] :
    const known_pos_range2 pr(9,20);

    BOOST_REQUIRE(! pr.is_pos_intersect(8));
    BOOST_REQUIRE(  pr.is_pos_intersect(9));
    BOOST_REQUIRE(  pr.is_pos_intersect(19));
    BOOST_REQUIRE(! pr.is_pos_intersect(20));
}


BOOST_AUTO_TEST_CASE( test_pos_range_semibound_is_pos_intersect )
{
    // this corresponds to zero-index range [-inf,19] :
    known_pos_range2 pr;
    pr.set_end_pos(20);

    BOOST_REQUIRE(  pr.is_pos_intersect(8));
    BOOST_REQUIRE(  pr.is_pos_intersect(9));
    BOOST_REQUIRE(  pr.is_pos_intersect(19));
    BOOST_REQUIRE(! pr.is_pos_intersect(20));
}

BOOST_AUTO_TEST_CASE( test_known_pos_range2_is_range_intersect )
{
    // this corresponds to zero-index range [9,19] :
    const known_pos_range2 pr(9,20);

    // left-side:
    BOOST_REQUIRE(! pr.is_range_intersect(known_pos_range2(0,9)));
    BOOST_REQUIRE(  pr.is_range_intersect(known_pos_range2(0,10)));

    // right side:
    BOOST_REQUIRE(  pr.is_range_intersect(known_pos_range2(19,30)));
    BOOST_REQUIRE(! pr.is_range_intersect(known_pos_range2(20,30)));

    // superset:
    BOOST_REQUIRE(  pr.is_range_intersect(known_pos_range2(0,30)));

    // subset:
    BOOST_REQUIRE(  pr.is_range_intersect(known_pos_range2(12,15)));
}


BOOST_AUTO_TEST_CASE( test_known_pos_range2_intersect_window )
{
    // this corresponds to zero-index range [9,19] :
    const known_pos_range2 pr1(9,20);
    const known_pos_range2 pr2(30,40);

    BOOST_REQUIRE(! is_intersect_window(pr1, pr2));
    BOOST_REQUIRE(! is_intersect_window(pr2, pr1));

    BOOST_REQUIRE(! is_intersect_window(pr1, pr2, 10));
    BOOST_REQUIRE(! is_intersect_window(pr2, pr1, 10));

    BOOST_REQUIRE(is_intersect_window(pr1, pr2, 11));
    BOOST_REQUIRE(is_intersect_window(pr2, pr1, 11));
}


BOOST_AUTO_TEST_SUITE_END()
