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

#include "boost/test/unit_test.hpp"

#include "htsapi/bam_header_info.hpp"
#include "test/testAlignmentDataUtil.hpp"

BOOST_AUTO_TEST_SUITE(bam_header_info_test_suite)

/// Create a bam_header_info object from an empty hts header
BOOST_AUTO_TEST_CASE(test_bam_header_info_empty)
{
  HtslibBamHeaderManager emptyHtsHeader;
  bam_header_info        emptyBamHeader(emptyHtsHeader.get());
  BOOST_REQUIRE(emptyBamHeader.empty());
}

/// Create a bam header_info object with one chromosome
BOOST_AUTO_TEST_CASE(test_bam_header_info_1chrom)
{
  std::vector<bam_header_info::chrom_info> chromData;
  chromData.emplace_back("chr1", 1000);

  HtslibBamHeaderManager oneHtsHeader(chromData);
  bam_header_info        oneBamHeader(oneHtsHeader.get());

  BOOST_REQUIRE_EQUAL(oneBamHeader.chrom_data.size(), 1);
  BOOST_REQUIRE_EQUAL(oneBamHeader.chrom_to_index.size(), 1);
}

/// Create a bam header_info object with two chromosomes
BOOST_AUTO_TEST_CASE(test_bam_header_info_2chrom)
{
  std::vector<bam_header_info::chrom_info> chromData;
  chromData.emplace_back("chr1", 1000);
  chromData.emplace_back("chr2", 1000);

  HtslibBamHeaderManager twoHtsHeader(chromData);
  bam_header_info        twoBamHeader(twoHtsHeader.get());

  BOOST_REQUIRE_EQUAL(twoBamHeader.chrom_data.size(), 2);
  BOOST_REQUIRE_EQUAL(twoBamHeader.chrom_to_index.size(), 2);
}

BOOST_AUTO_TEST_SUITE_END()
