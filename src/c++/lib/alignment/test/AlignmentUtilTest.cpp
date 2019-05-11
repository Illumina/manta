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
