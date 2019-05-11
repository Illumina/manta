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

#include "GlobalJumpAligner.hpp"

#include "blt_util/align_path.hpp"

#include <string>

BOOST_AUTO_TEST_SUITE(test_GlobalJumpAligner)

typedef short int score_t;

static JumpAlignmentResult<score_t> testAlign(
    const std::string& seq, const std::string& ref1, const std::string& ref2)
{
  AlignmentScores<score_t>     scores(2, -4, -5, -1, -1);
  score_t                      jumpScore(-3);
  GlobalJumpAligner<score_t>   aligner(scores, jumpScore);
  JumpAlignmentResult<score_t> result;
  aligner.align(seq.begin(), seq.end(), ref1.begin(), ref1.end(), ref2.begin(), ref2.end(), result);

  return result;
}

static JumpAlignmentResult<score_t> testAlign2(
    const std::string& seq, const std::string& ref1, const std::string& ref2)
{
  static const AlignmentScores<score_t> scores(2, -4, -10, -1, -1);
  static const int                      jumpScore(-20);
  GlobalJumpAligner<score_t>            aligner(scores, jumpScore);
  JumpAlignmentResult<score_t>          result;
  aligner.align(seq.begin(), seq.end(), ref1.begin(), ref1.end(), ref2.begin(), ref2.end(), result);

  return result;
}

static JumpAlignmentResult<score_t> testAlign3(
    const std::string& seq, const std::string& ref1, const std::string& ref2)
{
  static const AlignmentScores<score_t> scores(2, -4, -2, 0, -1);
  static const int                      jumpScore(-20);
  GlobalJumpAligner<score_t>            aligner(scores, jumpScore);
  JumpAlignmentResult<score_t>          result;
  aligner.align(seq.begin(), seq.end(), ref1.begin(), ref1.end(), ref2.begin(), ref2.end(), result);

  return result;
}

BOOST_AUTO_TEST_CASE(test_GlobalJumpAligner0)
{
  static const std::string seq("ABABACDCDC");
  static const std::string ref1("ABABA");
  static const std::string ref2("CDCDC");

  JumpAlignmentResult<score_t> result = testAlign(seq, ref1, ref2);

  BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align1.apath), "5=");
  BOOST_REQUIRE_EQUAL(result.align1.beginPos, 0);
  BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align2.apath), "5=");
  BOOST_REQUIRE_EQUAL(result.align2.beginPos, 0);
}

BOOST_AUTO_TEST_CASE(test_GlobalJumpAligner1)
{
  static const std::string seq("ABABACDCDC");
  static const std::string ref1("ABABAX");
  static const std::string ref2("CDCDC");

  JumpAlignmentResult<score_t> result = testAlign(seq, ref1, ref2);

  BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align1.apath), "5=");
  BOOST_REQUIRE_EQUAL(result.align1.beginPos, 0);
  BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align2.apath), "5=");
  BOOST_REQUIRE_EQUAL(result.align2.beginPos, 0);
}

BOOST_AUTO_TEST_CASE(test_GlobalJumpAligner2)
{
  static const std::string seq("ABABACDCDC");
  static const std::string ref1("ABABA");
  static const std::string ref2("XCDCDC");

  JumpAlignmentResult<score_t> result = testAlign(seq, ref1, ref2);

  BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align1.apath), "5=");
  BOOST_REQUIRE_EQUAL(result.align1.beginPos, 0);
  BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align2.apath), "5=");
  BOOST_REQUIRE_EQUAL(result.align2.beginPos, 1);
}

BOOST_AUTO_TEST_CASE(test_GlobalJumpAlignerLong)
{
  static const std::string seq("ABABACDCDC");
  static const std::string ref1("dslfjfkjaslABABAlsjfkdsflsk");
  static const std::string ref2("sdfldsklkjdCDCDCfsdlkjfslk");

  JumpAlignmentResult<score_t> result = testAlign(seq, ref1, ref2);

  BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align1.apath), "5=");
  BOOST_REQUIRE_EQUAL(result.align1.beginPos, 11);
  BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align2.apath), "5=");
  BOOST_REQUIRE_EQUAL(result.align2.beginPos, 11);
}

BOOST_AUTO_TEST_CASE(test_GlobalJumpAlignerSimpleIndels)
{
  static const std::string seq("ABABAABABACDCDCDyCDCDC");
  static const std::string ref1("xABABABABABAx");
  static const std::string ref2("xCDCDCDCDCDCDCx");

  JumpAlignmentResult<score_t> result = testAlign(seq, ref1, ref2);

  BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align1.apath), "5=1D5=");
  BOOST_REQUIRE_EQUAL(result.align1.beginPos, 1);
  BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align2.apath), "6=1I5=");
  BOOST_REQUIRE_EQUAL(result.align2.beginPos, 1);
}

BOOST_AUTO_TEST_CASE(test_GlobalJumpAlignerBPInsert)
{
  static const std::string seq("ABABABABABA1234CDCDCDCDCDC");
  static const std::string ref1("xABABABABABAx");
  static const std::string ref2("xCDCDCDCDCDCDCx");

  JumpAlignmentResult<score_t> result = testAlign(seq, ref1, ref2);

  BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align1.apath), "11=");
  BOOST_REQUIRE_EQUAL(result.align1.beginPos, 1);
  BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align2.apath), "11=");
  BOOST_REQUIRE_EQUAL(result.align2.beginPos, 1);
  BOOST_REQUIRE_EQUAL(result.jumpInsertSize, 4u);
}

// define behavior when the breakpoint solution is a 1d range
BOOST_AUTO_TEST_CASE(test_GlobalJumpAlignerBPRange)
{
  static const std::string seq("xyzxyzxyzABCABCABCxyzxyzxyz");
  static const std::string ref1("xyzxyzxyzxyzABCABCstustu");
  static const std::string ref2("stustuABCABCxyzxyzxyzxyz");

  JumpAlignmentResult<score_t> result = testAlign(seq, ref1, ref2);

  BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align1.apath), "12=");
  BOOST_REQUIRE_EQUAL(result.align1.beginPos, 3);
  BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align2.apath), "15=");
  BOOST_REQUIRE_EQUAL(result.align2.beginPos, 6);
  BOOST_REQUIRE_EQUAL(result.jumpInsertSize, 0u);
  BOOST_REQUIRE_EQUAL(result.jumpRange, 3u);
}

// define behavior when the breakpoint solution is a 1d range
BOOST_AUTO_TEST_CASE(test_GlobalJumpAlignerBPRange2)
{
  static const std::string seq("xyzxyzxyzABCABCABCxyzxyzxyz");
  static const std::string ref1("xyzxyzxyzxyzABCABCABCABCABC");
  static const std::string ref2("ABCABCABCABCABCxyzxyzxyzxyz");

  JumpAlignmentResult<score_t> result = testAlign(seq, ref1, ref2);

  BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align1.apath), "9=");
  BOOST_REQUIRE_EQUAL(result.align1.beginPos, 3);
  BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align2.apath), "18=");
  BOOST_REQUIRE_EQUAL(result.align2.beginPos, 6);
  BOOST_REQUIRE_EQUAL(result.jumpInsertSize, 0u);
  BOOST_REQUIRE_EQUAL(result.jumpRange, 9u);
}

// define behavior when the breakpoint solution has an insertion with a repeat
BOOST_AUTO_TEST_CASE(test_GlobalJumpAlignerInsert)
{
  static const std::string seq("xyzxyzxyzABCABCABCABCABCABCxyzxyzxyz");
  static const std::string ref1("xyzxyzxyzxyzABCABCstustu");
  static const std::string ref2("stustuABCABCxyzxyzxyzxyz");

  JumpAlignmentResult<score_t> result = testAlign(seq, ref1, ref2);

  BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align1.apath), "15=");
  BOOST_REQUIRE_EQUAL(result.align1.beginPos, 3);
  BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align2.apath), "15=");
  BOOST_REQUIRE_EQUAL(result.align2.beginPos, 6);
  BOOST_REQUIRE_EQUAL(result.jumpInsertSize, 6u);
  BOOST_REQUIRE_EQUAL(result.jumpRange, 0u);
}

BOOST_AUTO_TEST_CASE(test_GlobalJumpAlignerOnly1)
{
  static const std::string seq("ABABA");
  static const std::string ref1("xABABAx");
  static const std::string ref2("xCDCDCx");

  JumpAlignmentResult<score_t> result = testAlign(seq, ref1, ref2);

  BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align1.apath), "5=");
  BOOST_REQUIRE_EQUAL(result.align1.beginPos, 1);
  BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align2.apath), "");
  BOOST_REQUIRE_EQUAL(result.align2.beginPos, 0);
}

BOOST_AUTO_TEST_CASE(test_GlobalJumpAlignerOnly2)
{
  static const std::string seq("CDCDC");
  static const std::string ref1("xABABAx");
  static const std::string ref2("xCDCDCx");

  JumpAlignmentResult<score_t> result = testAlign(seq, ref1, ref2);

  BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align1.apath), "");
  BOOST_REQUIRE_EQUAL(result.align1.beginPos, 0);
  BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align2.apath), "5=");
  BOOST_REQUIRE_EQUAL(result.align2.beginPos, 1);
}

BOOST_AUTO_TEST_CASE(test_GlobalJumpAlignerOffEdge)
{
  static const std::string seq("123456ABABACDCDC123456");
  static const std::string ref1("xABABAx");
  static const std::string ref2("xCDCDCx");

  JumpAlignmentResult<score_t> result = testAlign(seq, ref1, ref2);

  BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align1.apath), "5S1X5=");
  BOOST_REQUIRE_EQUAL(result.align1.beginPos, 0);
  BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align2.apath), "5=1X5S");
  BOOST_REQUIRE_EQUAL(result.align2.beginPos, 1);
}

BOOST_AUTO_TEST_CASE(test_GlobalJumpAlignerRef2Clip)
{
  // extracted from production failure case:
  //
  static const std::string seq("GGCAGAAAAGGAAATA");
  static const std::string ref1("TAAAAAGTAGAT");
  static const std::string ref2("AAAGGAAATA");

  JumpAlignmentResult<score_t> result = testAlign2(seq, ref1, ref2);

  BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align1.apath), "");
  BOOST_REQUIRE_EQUAL(result.align1.beginPos, 0);
  BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align2.apath), "6S10=");
  BOOST_REQUIRE_EQUAL(result.align2.beginPos, 0);
  BOOST_REQUIRE_EQUAL(result.jumpInsertSize, 0u);
  BOOST_REQUIRE_EQUAL(result.jumpRange, 0u);
}

BOOST_AUTO_TEST_CASE(test_GlobalJumpAlignerRef1Clip)
{
  // extracted from production failure case:
  //
  static const std::string seq("TAAAAAGTAGATTTCGT");
  static const std::string ref1("TAAAAAGTAGAT");
  static const std::string ref2("AAAGGAAATA");

  JumpAlignmentResult<score_t> result = testAlign2(seq, ref1, ref2);

  BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align1.apath), "12=5S");
  BOOST_REQUIRE_EQUAL(result.align1.beginPos, 0);
  BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align2.apath), "");
  BOOST_REQUIRE_EQUAL(result.align2.beginPos, 0);
  BOOST_REQUIRE_EQUAL(result.jumpInsertSize, 0u);
  BOOST_REQUIRE_EQUAL(result.jumpRange, 0u);
}

BOOST_AUTO_TEST_CASE(test_GlobalJumpAlignerDelJunction)
{
  static const std::string seq("AAAACCCCCCCCTTTTAAAATTTT");
  static const std::string ref1("AAAAAAAACCCCCCCCG");
  static const std::string ref2("GGGGTTTTAAAATTTT");

  JumpAlignmentResult<score_t> result = testAlign3(seq, ref1, ref2);

  BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align1.apath), "12=");
  BOOST_REQUIRE_EQUAL(result.align1.beginPos, 4);
  BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align2.apath), "12=");
  BOOST_REQUIRE_EQUAL(result.align2.beginPos, 4);
  BOOST_REQUIRE_EQUAL(result.jumpInsertSize, 0u);
  BOOST_REQUIRE_EQUAL(result.jumpRange, 0u);
}

// define behavior of the breakpoint homology with and without a breakpoint insert
BOOST_AUTO_TEST_CASE(test_GlobalJumpAlignerJumpRange)
{
  static const std::string ref1("xxxxxxxxxABCABC");
  static const std::string ref2("ABCABCAyyyyyyyyy");

  // test1: the breakpoint has a homology seq "A", but no insert seq
  static const std::string     seq1("xxxxxxAyyyyyy");
  JumpAlignmentResult<score_t> result1 = testAlign(seq1, ref1, ref2);

  BOOST_REQUIRE_EQUAL(result1.align1.beginPos, 3);
  BOOST_REQUIRE_EQUAL(result1.align2.beginPos, 6);
  BOOST_REQUIRE_EQUAL(result1.jumpInsertSize, 0u);
  BOOST_REQUIRE_EQUAL(result1.jumpRange, 1u);

  // test2: the breakpoint has a insert seq "T", but no homology seq
  static const std::string     seq2("xxxxxxTAyyyyyy");
  JumpAlignmentResult<score_t> result2 = testAlign(seq2, ref1, ref2);

  BOOST_REQUIRE_EQUAL(result2.align1.beginPos, 3);
  BOOST_REQUIRE_EQUAL(result2.align2.beginPos, 6);
  BOOST_REQUIRE_EQUAL(result2.jumpInsertSize, 1u);
  BOOST_REQUIRE_EQUAL(result2.jumpRange, 0u);

  // test3: the breakpoint has a insert seq "T", but no homology seq
  static const std::string     seq3("xxxxxxATyyyyyy");
  JumpAlignmentResult<score_t> result3 = testAlign(seq3, ref1, ref2);

  BOOST_REQUIRE_EQUAL(result3.align1.beginPos, 3);
  BOOST_REQUIRE_EQUAL(result3.align2.beginPos, 7);
  BOOST_REQUIRE_EQUAL(result3.jumpInsertSize, 1u);
  BOOST_REQUIRE_EQUAL(result3.jumpRange, 0u);
}

BOOST_AUTO_TEST_SUITE_END()
