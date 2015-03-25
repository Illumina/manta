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

#include "blt_util/SimpleAlignment.hh"


BOOST_AUTO_TEST_SUITE( SimpleAlignment_path )


BOOST_AUTO_TEST_CASE( test_SimpleAlignment )
{
    using namespace ALIGNPATH;

    SimpleAlignment al;

    al.pos = 100;
    cigar_to_apath("100M", al.path);

    const known_pos_range2 testRange(matchifyEdgeSoftClipRefRange(al));

    BOOST_REQUIRE_EQUAL(testRange.begin_pos(), 100);
    BOOST_REQUIRE_EQUAL(testRange.end_pos(), 200);
}


BOOST_AUTO_TEST_CASE( test_SimpleAlignment2 )
{
    using namespace ALIGNPATH;

    SimpleAlignment al;

    al.pos = 100;
    cigar_to_apath("10S50M10D40M10S", al.path);

    const known_pos_range2 testRange(matchifyEdgeSoftClipRefRange(al));

    BOOST_REQUIRE_EQUAL(testRange.begin_pos(), 90);
    BOOST_REQUIRE_EQUAL(testRange.end_pos(), 210);
}


BOOST_AUTO_TEST_SUITE_END()

