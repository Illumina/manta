//
// Manta - Structural Variant and Indel Caller
// Copyright (c) 2013-2018 Illumina, Inc.
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

#include "testConfig.h"

#include "htsapi/bed_streamer.hh"

#include "boost/test/unit_test.hpp"



BOOST_AUTO_TEST_SUITE( test_bed_streamer )

static
const char*
getTestpath()
{
    static const std::string testPath(std::string(TEST_DATA_PATH) + "/bed_streamer_test.bed.gz");
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

