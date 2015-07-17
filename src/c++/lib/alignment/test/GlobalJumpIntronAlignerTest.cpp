// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Manta - Structural Variant and Indel Caller
// Copyright (c) 2013-2015 Illumina, Inc.
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

///
/// \author Felix Schlesinger
///

#include "boost/test/unit_test.hpp"

#include "GlobalJumpIntronAligner.hh"

#include "blt_util/align_path.hh"

#include <string>



BOOST_AUTO_TEST_SUITE( test_GlobalJumpIntronAligner )

typedef short int score_t;

static
JumpAlignmentResult<score_t>
testAlignScores(
    const std::string& seq,
    const std::string& ref1,
    const std::string& ref2,
    int match, int mismatch, int open, int extend, int spliceOpen, int offEdge, int spliceOffEdge,
    bool stranded, bool bp1Fw, bool bp2Fw)
{
    AlignmentScores<score_t> scores(match, mismatch, open, extend, offEdge);
    score_t jumpScore(-3);
    GlobalJumpIntronAligner<score_t> aligner(scores,jumpScore,spliceOpen,spliceOffEdge);
    JumpAlignmentResult<score_t> result;
    aligner.align(
        seq.begin(),seq.end(),
        ref1.begin(),ref1.end(),
        ref2.begin(),ref2.end(),
        bp1Fw, bp2Fw, stranded,
        result);

    return result;
}

static
JumpAlignmentResult<score_t>
testAlignSplice(
    const std::string& seq,
    const std::string& ref1,
    const std::string& ref2,
    bool stranded=true,
    bool fw=true)
{
    return testAlignScores(seq, ref1, ref2, 2,-4,-5,-1,-4,-1,-1, stranded, fw, fw);
}


BOOST_AUTO_TEST_CASE( test_GlobalJumpAlignerSplice )
{
    static const std::string seq("AAAAABBBBB");
    static const std::string ref1("xAAAAAGTxxxAGBBBBBx");
    static const std::string ref2("xxxx");

    JumpAlignmentResult<score_t> result = testAlignSplice(seq,ref1,ref2);

    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align1.apath),"5=7N5=");
    BOOST_REQUIRE_EQUAL(result.align1.beginPos,1);
}

BOOST_AUTO_TEST_CASE( test_GlobalJumpAlignerSpliceRef2 )
{
    static const std::string seq("AAAAABBBBB");
    static const std::string ref1("xxxx");
    static const std::string ref2("xAAAAAGTxxxAGBBBBBx");

    JumpAlignmentResult<score_t> result = testAlignSplice(seq,ref1,ref2);

    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align2.apath),"5=7N5=");
    BOOST_REQUIRE_EQUAL(result.align2.beginPos,1);
}

BOOST_AUTO_TEST_CASE(test_GlobalJumpAlignerSpliceRev)
{
    static const std::string seq("AAAAABBBBB");
    static const std::string ref1("xxxx");
    static const std::string ref2("xAAAAACTxxxACBBBBBx");

    JumpAlignmentResult<score_t> result = testAlignSplice(seq, ref1, ref2, true, false);

    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align2.apath), "5=7N5=");
    BOOST_REQUIRE_EQUAL(result.align2.beginPos, 1);
}

BOOST_AUTO_TEST_CASE(test_GlobalJumpAlignerSpliceUnstrandedRev)
{
    static const std::string seq("AAAAABBBBB");
    static const std::string ref1("xxxx");
    static const std::string ref2("xAAAAACTxxxACBBBBBx");

    JumpAlignmentResult<score_t> result = testAlignSplice(seq, ref1, ref2, false);

    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align2.apath), "5=7N5=");
    BOOST_REQUIRE_EQUAL(result.align2.beginPos, 1);
}

BOOST_AUTO_TEST_CASE(test_GlobalJumpAlignerNoSpliceWrongStrand)
{
    static const std::string seq("AAAAABBBBB");
    static const std::string ref1("xxxx");
    static const std::string ref2("xAAAAACTxxxACBBBBBx");

    JumpAlignmentResult<score_t> result = testAlignSplice(seq, ref1, ref2, true, true);

    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align2.apath), "5=7D5=");
    BOOST_REQUIRE_EQUAL(result.align2.beginPos, 1);
}

BOOST_AUTO_TEST_CASE( test_GlobalJumpAlignerSpliceNoSplice )
{
    static const std::string seq("AAAAABBBBB");
    static const std::string ref1("xAAAAAGGxxxAGBBBBBx");
    static const std::string ref2("xxxx");

    JumpAlignmentResult<score_t> result = testAlignSplice(seq,ref1,ref2);

    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align1.apath),"5=7D5=");
    BOOST_REQUIRE_EQUAL(result.align1.beginPos,1);
}

BOOST_AUTO_TEST_CASE( test_GlobalJumpAlignerSpliceNoSplice2 )
{
    static const std::string seq("AAAAABBBBB");
    static const std::string ref1("xAAAAAGTxxxGGBBBBBx");
    static const std::string ref2("xxxx");

    JumpAlignmentResult<score_t> result = testAlignSplice(seq,ref1,ref2);

    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align1.apath),"5=7D5=");
    BOOST_REQUIRE_EQUAL(result.align1.beginPos,1);
}

BOOST_AUTO_TEST_SUITE_END()
