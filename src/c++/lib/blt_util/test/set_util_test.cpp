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

#include "blt_util/set_util.hh"


BOOST_AUTO_TEST_SUITE( test_set_util )


BOOST_AUTO_TEST_CASE( test_set_subtract )
{
    std::set<int> A;
    A.insert(1);
    A.insert(4);
    A.insert(5);
    A.insert(6);
    A.insert(9);

    std::set<int> B;
    B.insert(2);
    B.insert(5);
    B.insert(6);
    B.insert(7);
    B.insert(10);

    inplaceSetSubtract(A,B);

    BOOST_REQUIRE_EQUAL(B.size(), 3);
    BOOST_REQUIRE_EQUAL(*B.begin(), 2);
}


BOOST_AUTO_TEST_SUITE_END()

