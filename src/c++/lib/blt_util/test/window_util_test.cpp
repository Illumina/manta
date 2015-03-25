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

#include "window_util.hh"

#include <iostream>


BOOST_AUTO_TEST_SUITE( test_window )


BOOST_AUTO_TEST_CASE( test_window )
{
    static const double tol(0.0001);

    window_average wa(3);
    wa.insert(0);
    BOOST_REQUIRE_EQUAL(wa.size(),1);
    BOOST_REQUIRE_CLOSE(wa.avg(), 0., tol);
    wa.insert(1);
    BOOST_REQUIRE_EQUAL(wa.size(),2);
    BOOST_REQUIRE_CLOSE(wa.avg(), 0.5, tol);
    wa.insert(2);
    BOOST_REQUIRE_EQUAL(wa.size(),3);
    BOOST_REQUIRE_CLOSE(wa.avg(), 1., tol);
    wa.insert(3);
    BOOST_REQUIRE_EQUAL(wa.size(),3);
    BOOST_REQUIRE_CLOSE(wa.avg(), 2., tol);

    wa.insert_null();
    BOOST_REQUIRE_EQUAL(wa.size(),2);
    BOOST_REQUIRE_CLOSE(wa.avg(), 2.5, tol);
    wa.insert_null();
    BOOST_REQUIRE_EQUAL(wa.size(),1);
    BOOST_REQUIRE_CLOSE(wa.avg(), 3., tol);
    wa.insert_null();
    BOOST_REQUIRE_EQUAL(wa.size(),0);

    wa.insert(0);
    BOOST_REQUIRE_EQUAL(wa.size(),1);
    BOOST_REQUIRE_CLOSE(wa.avg(), 0., tol);
    wa.insert(1);
    BOOST_REQUIRE_EQUAL(wa.size(),2);
    BOOST_REQUIRE_CLOSE(wa.avg(), 0.5, tol);
    wa.insert(2);
    BOOST_REQUIRE_EQUAL(wa.size(),3);
    BOOST_REQUIRE_CLOSE(wa.avg(), 1., tol);
    wa.insert(3);
    BOOST_REQUIRE_EQUAL(wa.size(),3);
    BOOST_REQUIRE_CLOSE(wa.avg(), 2., tol);
}

BOOST_AUTO_TEST_SUITE_END()

