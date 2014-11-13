// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Starka
// Copyright (c) 2009-2014 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/sequencing/licenses/>
//

#include "htsapi/bed_streamer.hh"

#include "boost/test/unit_test.hpp"



BOOST_AUTO_TEST_SUITE( test_bed_streamer )

static
const char*
getTestpath()
{
    static const std::string testPath("@CMAKE_CURRENT_SOURCE_DIR@/bed_streamer_test.bed.gz");
    return testPath.c_str();
}


BOOST_AUTO_TEST_CASE( test_bed_streamer_region )
{
    bed_streamer beds(getTestpath(),"chr1:750000-822000");

    const bed_record* bedr(nullptr);

    BOOST_REQUIRE( beds.next() );
    bedr = beds.get_record_ptr();
    assert(bedr != nullptr);

    BOOST_REQUIRE( bedr->is_valid() );

    BOOST_REQUIRE_EQUAL(bedr->chrom, "chr1");
    BOOST_REQUIRE_EQUAL(bedr->begin, 750000);
    BOOST_REQUIRE_EQUAL(bedr->end, 750001);

    BOOST_REQUIRE( beds.next() );
    bedr = beds.get_record_ptr();
    assert(bedr != nullptr);

    BOOST_REQUIRE( bedr->is_valid() );
    BOOST_REQUIRE_EQUAL(bedr->chrom, "chr1");
    BOOST_REQUIRE_EQUAL(bedr->begin, 800000);
    BOOST_REQUIRE_EQUAL(bedr->end, 800001);

    BOOST_REQUIRE(! beds.next() );
}


BOOST_AUTO_TEST_SUITE_END()

