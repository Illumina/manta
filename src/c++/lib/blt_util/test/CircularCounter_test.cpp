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

#include "CircularCounter.hh"


BOOST_AUTO_TEST_SUITE( test_CircularCounter )


BOOST_AUTO_TEST_CASE( test_CircularCounter1 )
{
    CircularCounter cc(3);

    BOOST_CHECK_EQUAL(cc.count(),0);

    cc.push(true);
    BOOST_CHECK_EQUAL(cc.count(),1);

    cc.push(false);
    cc.push(false);
    BOOST_CHECK_EQUAL(cc.count(),1);
    cc.push(false);
    BOOST_CHECK_EQUAL(cc.count(),0);
}


BOOST_AUTO_TEST_SUITE_END()
