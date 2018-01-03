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

#include "test_config.h"

#include "htsapi/bam_streamer.hh"

#include "boost/test/unit_test.hpp"



BOOST_AUTO_TEST_SUITE( test_bam_streamer )


BOOST_AUTO_TEST_CASE( test_bam_streamer_bam_read )
{
    const std::string testBamPath(std::string(TEST_DATA_PATH) + "/alignment_test.bam");

    /// Assert that reference pointer is optional for BAM
    bam_streamer stream(testBamPath.c_str(), nullptr);

    {
        // iterate through entire file:
        unsigned count(0);
        while (stream.next())
        {
            const bam_record& read(*(stream.get_record_ptr()));
            if (! read.is_unmapped()) count++;
        }
        BOOST_REQUIRE_EQUAL(count, 4u);
    }

    {
        // iterate through region:
        stream.resetRegion("chrA");
        unsigned count(0);
        while (stream.next())
        {
            const bam_record& read(*(stream.get_record_ptr()));
            if (! read.is_unmapped()) count++;
        }
        BOOST_REQUIRE_EQUAL(count, 2u);
    }
}


BOOST_AUTO_TEST_CASE( test_bam_streamer_cram_read )
{
    const std::string testCramPath(std::string(TEST_DATA_PATH) + "/alignment_test.cram");
    const std::string testRefPath(std::string(TEST_DATA_PATH) + "/alignment_test.fasta");

    bam_streamer stream(testCramPath.c_str(), testRefPath.c_str());

    {
        // iterate through entire file:
        unsigned count(0);
        while (stream.next())
        {
            const bam_record& read(*(stream.get_record_ptr()));
            if (! read.is_unmapped()) count++;
        }
        BOOST_REQUIRE_EQUAL(count, 4u);
    }

    {
        // iterate through region:
        stream.resetRegion("chrA");
        unsigned count(0);
        while (stream.next())
        {
            const bam_record& read(*(stream.get_record_ptr()));
            if (! read.is_unmapped()) count++;
        }
        BOOST_REQUIRE_EQUAL(count, 2u);
    }
}


BOOST_AUTO_TEST_SUITE_END()
