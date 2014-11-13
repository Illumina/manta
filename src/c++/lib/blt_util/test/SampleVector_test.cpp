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

#include "SampleVector.hh"

#include <random>


BOOST_AUTO_TEST_SUITE( test_SampleVector )


BOOST_AUTO_TEST_CASE( test_SampleVector1 )
{
    std::mt19937 rngEngine(0);
    SampleVector<int,std::mt19937> sv(2,rngEngine);

    for (unsigned i(0); i<100; ++i)
    {
        sv.push(i);
    }

    BOOST_CHECK_EQUAL(sv.data()[0],34);
    BOOST_CHECK_EQUAL(sv.data()[1],70);
}


BOOST_AUTO_TEST_SUITE_END()
