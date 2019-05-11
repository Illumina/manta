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

#include "htsapi/align_path_bam_util.hpp"
#include "htsapi/bam_record.hpp"

#include "boost/test/unit_test.hpp"

BOOST_AUTO_TEST_SUITE(test_align_path_bam_util)

BOOST_AUTO_TEST_CASE(test_edit_bam_cigar)
{
  const std::string testCigar("10M10D1I10M");
  ALIGNPATH::path_t inputPath;
  cigar_to_apath(testCigar.c_str(), inputPath);

  bam_record bamRead;
  bam1_t*    bamDataPtr(bamRead.get_data());
  edit_bam_cigar(inputPath, *bamDataPtr);

  ALIGNPATH::path_t outputPath;
  bam_cigar_to_apath(bamRead.raw_cigar(), bamRead.n_cigar(), outputPath);

  BOOST_REQUIRE_EQUAL(apath_to_cigar(outputPath), testCigar);
}

BOOST_AUTO_TEST_SUITE_END()
