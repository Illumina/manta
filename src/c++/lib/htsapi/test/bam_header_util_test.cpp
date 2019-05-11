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
