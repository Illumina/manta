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

#include "GlobalLargeIndelAligner.hpp"

#include "blt_util/align_path.hpp"

#include <string>

BOOST_AUTO_TEST_SUITE(test_GlobalLargeIndelAligner)

typedef short int score_t;

static AlignmentResult<score_t> testAlign(const std::string& seq, const std::string& ref)
{
  AlignmentScores<score_t>         scores(2, -4, -5, -1, -4);
  GlobalLargeIndelAligner<score_t> aligner(scores, -10);
  AlignmentResult<score_t>         result;
  aligner.align(seq.begin(), seq.end(), ref.begin(), ref.end(), result);

  return result;
}

BOOST_AUTO_TEST_CASE(test_GlobalLargeIndelAligner1)
{
  static const std::string seq("D");
  static const std::string ref("ABCDEF");

  AlignmentResult<score_t> result = testAlign(seq, ref);

  BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align.apath), "1=");
  BOOST_REQUIRE_EQUAL(result.align.beginPos, 3);
}

BOOST_AUTO_TEST_CASE(test_GlobalLargeIndelAlignerDelete)
{
  static const std::string seq("BCDEFHIKLM");
  static const std::string ref("ABCDEFGHIKLMN");

  AlignmentResult<score_t> result = testAlign(seq, ref);

  BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align.apath), "5=1D5=");
  BOOST_REQUIRE_EQUAL(result.align.beginPos, 1);
}

BOOST_AUTO_TEST_CASE(test_GlobalLargeIndelAlignerInsert)
{
  static const std::string seq("BCDEFGXHIKLM");
  static const std::string ref("ABCDEFGHIKLMN");

  AlignmentResult<score_t> result = testAlign(seq, ref);

  BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align.apath), "6=1I5=");
  BOOST_REQUIRE_EQUAL(result.align.beginPos, 1);
}

BOOST_AUTO_TEST_CASE(test_GlobalLargeIndelAlignerInsertDelete)
{
  static const std::string seq("BBBBBBCDXYZHIKLMMMM");
  static const std::string ref("ABBBBBBCDEFGHIKLMMMMN");

  AlignmentResult<score_t> result = testAlign(seq, ref);

  BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align.apath), "8=3I3D8=");
  BOOST_REQUIRE_EQUAL(result.align.beginPos, 1);
}

BOOST_AUTO_TEST_CASE(test_GlobalLargeIndelAlignerInsertDelete2)
{
  static const std::string seq("BBBBBBCDEXYHIKLMMMM");
  static const std::string ref("ABBBBBBCDEFGHIKLMMMMN");

  AlignmentResult<score_t> result = testAlign(seq, ref);

  BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align.apath), "9=2X8=");
  BOOST_REQUIRE_EQUAL(result.align.beginPos, 1);
}

BOOST_AUTO_TEST_CASE(test_GlobalLargeIndelAlignerShortRef1)
{
  static const std::string seq("ABCD");
  static const std::string ref("BCD");

  AlignmentResult<score_t> result = testAlign(seq, ref);

  BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align.apath), "1S3=");
  BOOST_REQUIRE_EQUAL(result.align.beginPos, 0);
  BOOST_REQUIRE_EQUAL(result.score, 2);
}

BOOST_AUTO_TEST_CASE(test_GlobalLargeIndelAlignerShortRef2)
{
  static const std::string seq("ABCD");
  static const std::string ref("ABC");

  AlignmentResult<score_t> result = testAlign(seq, ref);

  BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align.apath), "3=1S");
  BOOST_REQUIRE_EQUAL(result.align.beginPos, 0);
  BOOST_REQUIRE_EQUAL(result.score, 2);
}

BOOST_AUTO_TEST_CASE(test_GlobalLargeIndelAlignerShortRef3)
{
  static const std::string seq("ABCD");
  static const std::string ref("B");

  AlignmentResult<score_t> result = testAlign(seq, ref);

  BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align.apath), "1S1=2S");
  BOOST_REQUIRE_EQUAL(result.align.beginPos, 0);
  BOOST_REQUIRE_EQUAL(result.score, -10);
}

// show that the method left aligns a deletion within a repeat
//
BOOST_AUTO_TEST_CASE(test_GlobalLargeIndelAlignerLeftShift)
{
  static const std::string seq("ABCDEFFFFFGHIJKL");
  static const std::string ref("ABCDEFFFFFFGHIJKL");

  AlignmentResult<score_t> result = testAlign(seq, ref);

  BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align.apath), "5=1D11=");
  BOOST_REQUIRE_EQUAL(result.align.beginPos, 0);
}

// show that the method left aligns an insertion within a repeat
//
BOOST_AUTO_TEST_CASE(test_GlobalLargeIndelAlignerLeftShift2)
{
  static const std::string seq("ABCDEFFFFFFFGHIJKL");
  static const std::string ref("ABCDEFFFFFFGHIJKL");

  AlignmentResult<score_t> result = testAlign(seq, ref);

  BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align.apath), "5=1I12=");
  BOOST_REQUIRE_EQUAL(result.align.beginPos, 0);
}

BOOST_AUTO_TEST_CASE(test_GlobalLargeIndelAlignerBigDelete)
{
  static const std::string seq("BCDEFHIKLM");
  static const std::string ref("ABCDEFGGGGGGGGGGGGGGGGGGGGGGGGGGHIKLMN");

  AlignmentResult<score_t> result = testAlign(seq, ref);

  BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align.apath), "5=26D5=");
  BOOST_REQUIRE_EQUAL(result.align.beginPos, 1);
  BOOST_REQUIRE_EQUAL(result.score, 10);
}

BOOST_AUTO_TEST_CASE(test_GlobalLargeIndelAlignerBigDeleteSmallInsert)
{
  static const std::string seq("BCDEFXHIKLM");
  static const std::string ref("ABCDEFGGGGGGGGGGGGGGGGGGGGGGGGGGHIKLMN");

  AlignmentResult<score_t> result = testAlign(seq, ref);

  BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align.apath), "5=1I26D5=");
  BOOST_REQUIRE_EQUAL(result.align.beginPos, 1);
  BOOST_REQUIRE_EQUAL(result.score, 9);
}

BOOST_AUTO_TEST_CASE(test_GlobalLargeIndelAlignerBigInsert)
{
  static const std::string seq("BCDEFGXXXXXXXXXXXXXXXXXXXXXXXXHIKLM");
  static const std::string ref("ABCDEFGHIKLMN");

  AlignmentResult<score_t> result = testAlign(seq, ref);

  BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align.apath), "6=24I5=");
  BOOST_REQUIRE_EQUAL(result.score, 12);
}

BOOST_AUTO_TEST_CASE(test_GlobalLargeIndelAlignerBigInsertSmallDelete)
{
  static const std::string seq("BCDEFGXXXXXXXXXXXXXXXXXXXXXXXXIKLM");
  static const std::string ref("ABCDEFGHIKLMN");

  AlignmentResult<score_t> result = testAlign(seq, ref);

  BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align.apath), "6=24I1D4=");
  BOOST_REQUIRE_EQUAL(result.score, 9);
}

BOOST_AUTO_TEST_CASE(test_GlobalLargeIndelAlignerBigInsertSmallDelete2)
{
  static const std::string seq("BCDEFXXXXXXXXXXXXXXXXXXXXXXXXHIKLM");
  static const std::string ref("ABCDEFGHIKLMN");

  AlignmentResult<score_t> result = testAlign(seq, ref);

  BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align.apath), "5=24I1D5=");
  BOOST_REQUIRE_EQUAL(result.score, 9);
}

BOOST_AUTO_TEST_SUITE_END()
