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
