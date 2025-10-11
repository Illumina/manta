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

#include "common/Exceptions.hpp"
#include "htsapi/bed_streamer.hpp"

#include "boost/test/unit_test.hpp"

static const char* getTestpath()
{
  static const std::string testPath(std::string(TEST_DATA_PATH) + "/bed_streamer_test.bed.gz");
  return testPath.c_str();
}

static const char* getTestpath2()
{
  static const std::string testPath(std::string(TEST_DATA_PATH) + "/bed_streamer_test2.bed.gz");
  return testPath.c_str();
}

BOOST_AUTO_TEST_SUITE(bed_streamer_test_suite)

BOOST_AUTO_TEST_CASE(test_bed_streamer_region)
{
  bed_streamer beds(getTestpath(), "chr1:750000-822000");

  const bed_record* bedr(nullptr);

  BOOST_REQUIRE(beds.next());
  bedr = beds.get_record_ptr();
  assert(bedr != nullptr);

  BOOST_REQUIRE(bedr->is_valid());

  BOOST_REQUIRE_EQUAL(bedr->chrom, "chr1");
  BOOST_REQUIRE_EQUAL(bedr->begin, 750000);
  BOOST_REQUIRE_EQUAL(bedr->end, 750001);

  BOOST_REQUIRE(beds.next());
  bedr = beds.get_record_ptr();
  assert(bedr != nullptr);

  BOOST_REQUIRE(bedr->is_valid());
  BOOST_REQUIRE_EQUAL(bedr->chrom, "chr1");
  BOOST_REQUIRE_EQUAL(bedr->begin, 800000);
  BOOST_REQUIRE_EQUAL(bedr->end, 800001);

  BOOST_REQUIRE(!beds.next());
}

// Test behavior on an invalid bed record
BOOST_AUTO_TEST_CASE(test_bed_streamer_region2)
{
  {
    const bool   requireNonZeroRegionLength(true);
    bed_streamer beds(getTestpath2(), "demo20:1-5000", requireNonZeroRegionLength);
    BOOST_REQUIRE_THROW(beds.next(), illumina::common::GeneralException);
  }

  {
    const bool   requireNonZeroRegionLength(false);
    bed_streamer beds(getTestpath2(), "demo20:1-5000", requireNonZeroRegionLength);
    BOOST_REQUIRE(!beds.next());
  }
}

BOOST_AUTO_TEST_SUITE_END()
