//
// Manta - Structural Variant and Indel Caller
// Copyright (c) 2013-2025 Illumina, Inc.
//
// This program is licensed under the terms of the Polyform strict license
//
// ***As far as the law allows, the software comes as is, without
// any warranty or condition, and the licensor will not be liable
// to you for any damages arising out of these terms or the use
// or nature of the software, under any kind of legal claim.***
//
// You should have received a copy of the PolyForm Strict License 1.0.0
// along with this program.  If not, see <https://polyformproject.org/licenses/strict/1.0.0>.
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
