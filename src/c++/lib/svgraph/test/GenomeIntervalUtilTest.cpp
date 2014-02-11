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

///
/// \author Chris Saunders
///

#include "boost/test/unit_test.hpp"

#include "svgraph/GenomeIntervalUtil.hh"



BOOST_AUTO_TEST_SUITE( test_GenomeIntervalUtil )

BOOST_AUTO_TEST_CASE( test_IntervalCompressor )
{
    // test that GenomeInterval sorting follows expect:
    std::vector<GenomeInterval> test;

    test.push_back(GenomeInterval(1,15,19));
    test.push_back(GenomeInterval(2,5,10));
    test.push_back(GenomeInterval(1,10,20));
    test.push_back(GenomeInterval(1,24,50));
    test.push_back(GenomeInterval(2,8,10));
    test.push_back(GenomeInterval(1,15,22));

    const std::vector<unsigned> indexMap = intervalCompressor(test);

    BOOST_REQUIRE_EQUAL(test.size(),3u);
    BOOST_REQUIRE_EQUAL(test[0],GenomeInterval(1,10,22));
    BOOST_REQUIRE_EQUAL(test[1],GenomeInterval(2,5,10));

    BOOST_REQUIRE_EQUAL(indexMap.size(),6u);
    BOOST_REQUIRE_EQUAL(indexMap[0],0u);
    BOOST_REQUIRE_EQUAL(indexMap[5],0u);
    BOOST_REQUIRE_EQUAL(indexMap[4],1u);
}


BOOST_AUTO_TEST_SUITE_END()
