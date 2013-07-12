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

///
/// \author Chris Saunders
///

#include "boost/test/unit_test.hpp"

#include "svgraph/GenomeInterval.hh"



BOOST_AUTO_TEST_SUITE( test_GenomeInterval )

BOOST_AUTO_TEST_CASE( test_GenomeInterval )
{

    // test that GenomeInterval sorting follows expect:
    std::vector<GenomeInterval> test;

    test.push_back(GenomeInterval(1,15,19));
    test.push_back(GenomeInterval(1,15,22));
    test.push_back(GenomeInterval(1,10,20));
    test.push_back(GenomeInterval(2,5,10));
    test.push_back(GenomeInterval(2,8,10));

    std::sort(test.begin(),test.end());

    BOOST_REQUIRE_EQUAL(test[0],GenomeInterval(1,10,20));
    BOOST_REQUIRE_EQUAL(test[2],GenomeInterval(1,15,22));
    BOOST_REQUIRE_EQUAL(test[4],GenomeInterval(2,8,10));
}

BOOST_AUTO_TEST_SUITE_END()

