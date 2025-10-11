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

#include "boost/test/unit_test.hpp"

#include "AlignmentUtil.hpp"

#include "blt_util/align_path.hpp"

#include <string>

BOOST_AUTO_TEST_SUITE(test_AlignmnetUtil)

BOOST_AUTO_TEST_CASE(test_getQuerySegments)
{
  const std::string querySeq("AAAACCGGG");

  JumpAlignmentResult<int> result;

  cigar_to_apath("4M", result.align1.apath);
  result.jumpInsertSize = 2;
  cigar_to_apath("3M", result.align2.apath);

  std::string bp1Seq, insertSeq, bp2Seq;
  getFwdStrandQuerySegments(result, querySeq, false, true, false, bp1Seq, bp2Seq, insertSeq);

  BOOST_REQUIRE_EQUAL(bp1Seq, "TTTT");
  BOOST_REQUIRE_EQUAL(insertSeq, "GG");
  BOOST_REQUIRE_EQUAL(bp2Seq, "GGG");
}

BOOST_AUTO_TEST_SUITE_END()
