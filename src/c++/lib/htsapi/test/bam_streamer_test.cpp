//
// Manta - Structural Variant and Indel Caller
// Copyright (c) 2013-2019 Illumina, Inc.
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

#include "blt_util/blt_exception.hpp"
#include "htsapi/bam_streamer.hpp"

#include "boost/test/unit_test.hpp"

BOOST_AUTO_TEST_SUITE(bam_streamer_test_suite)

static void checkStream(bam_streamer& stream, const unsigned expectedCount)
{
  unsigned count(0);
  while (stream.next()) {
    const bam_record& read(*(stream.get_record_ptr()));
    if (!read.is_unmapped()) count++;
  }
  BOOST_REQUIRE_EQUAL(count, expectedCount);
}

BOOST_AUTO_TEST_CASE(test_bam_streamer_bam_read)
{
  const std::string testBamPath(std::string(TEST_DATA_PATH) + "/alignment_test.bam");

  // Assert that reference pointer is optional for BAM
  bam_streamer stream(testBamPath.c_str(), nullptr);

  // iterate through entire file:
  checkStream(stream, 4u);

  // iterate through region:
  stream.resetRegion("chrA");
  checkStream(stream, 2u);
}

BOOST_AUTO_TEST_CASE(test_bam_streamer_cram_read)
{
  const std::string testCramPath(std::string(TEST_DATA_PATH) + "/alignment_test.cram");
  const std::string testRefPath(std::string(TEST_DATA_PATH) + "/alignment_test.fasta");

  bam_streamer stream(testCramPath.c_str(), testRefPath.c_str());

  // iterate through entire file:
  checkStream(stream, 4u);

  // iterate through region:
  stream.resetRegion("chrA");
  checkStream(stream, 2u);
}

BOOST_AUTO_TEST_CASE(test_bam_streamer_cram_read_fail)
{
  const std::string testCramPath(std::string(TEST_DATA_PATH) + "/alignment_test.cram");

  {
    bam_streamer stream(testCramPath.c_str(), nullptr);
    BOOST_REQUIRE_THROW(stream.next(), blt_exception);
  }

  {
    bam_streamer stream(testCramPath.c_str(), nullptr);
    stream.resetRegion("chrA");
    BOOST_REQUIRE_THROW(stream.next(), blt_exception);
  }
}

BOOST_AUTO_TEST_SUITE_END()
