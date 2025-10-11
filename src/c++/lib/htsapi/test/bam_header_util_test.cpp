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

#include "htsapi/bam_header_util.hpp"

#include "boost/test/unit_test.hpp"

BOOST_AUTO_TEST_SUITE(test_bam_header_util)

BOOST_AUTO_TEST_CASE(test_parse_bam_region)
{
  std::string chrom;
  int32_t     begin, end;
  parse_bam_region("HLA-A*01:01:01:02N:1-3291", chrom, begin, end);

  BOOST_REQUIRE_EQUAL(chrom, "HLA-A*01:01:01:02N");
  BOOST_REQUIRE_EQUAL(begin, 0);
  BOOST_REQUIRE_EQUAL(end, 3291);

  parse_bam_region("HLA-A*01:01:01:02N", chrom, begin, end);

  BOOST_REQUIRE_EQUAL(chrom, "HLA-A*01:01:01:02N");
  BOOST_REQUIRE_EQUAL(begin, 0);
  BOOST_REQUIRE(end > 1000000000);
}

BOOST_AUTO_TEST_SUITE_END()
