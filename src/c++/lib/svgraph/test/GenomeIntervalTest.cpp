//
// Manta - Structural Variant and Indel Caller
// Copyright (c) 2013-2017 Illumina, Inc.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//

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

