// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Manta
// Copyright (c) 2013-2014 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/sequencing/licenses/>
//

///
/// \author Chris Saunders
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
    int match, int mismatch, int open, int extend, int spliceOpen, int offEdge, int spliceOffEdge)
{
    AlignmentScores<score_t> scores(match, mismatch, open, extend, spliceOpen, offEdge, spliceOffEdge);
    score_t jumpScore(-3);
    GlobalJumpIntronAligner<score_t> aligner(scores,jumpScore);
    JumpAlignmentResult<score_t> result;
    aligner.align(
        seq.begin(),seq.end(),
        ref1.begin(),ref1.end(),
        ref2.begin(),ref2.end(),
        result);

    return result;
}

static
JumpAlignmentResult<score_t>
testAlignSplice(
    const std::string& seq,
    const std::string& ref1,
    const std::string& ref2)
{
    return testAlignScores(seq, ref1, ref2, 2,-4,-5,-1,-4,-1,-1);
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
