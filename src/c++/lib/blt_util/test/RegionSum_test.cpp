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

#include "RegionSum.hh"

BOOST_AUTO_TEST_SUITE( test_RegionSum )


BOOST_AUTO_TEST_CASE( RegionSum_test )
{
    RegionSum<unsigned> rs;
    rs.add(known_pos_range2(3,7),1u);
    rs.add(known_pos_range2(4,5),2u);

    BOOST_REQUIRE_EQUAL(rs.maxVal(),3u);
}


BOOST_AUTO_TEST_SUITE_END()
