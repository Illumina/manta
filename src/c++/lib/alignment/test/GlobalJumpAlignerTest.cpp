// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Manta
// Copyright (c) 2013 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/downloads/sequencing/licenses/>.
//

///
/// \author Chris Saunders
///

#include "boost/test/unit_test.hpp"

#include "GlobalJumpAligner.hh"

#include "blt_util/align_path.hh"

#include <string>



BOOST_AUTO_TEST_SUITE( test_GlobalJumpAligner )

typedef short int score_t;

static
JumpAlignmentResult<score_t>
testAlign(
    const std::string& seq,
    const std::string& ref1,
    const std::string& ref2)
{
    AlignmentScores<score_t> scores(2,-4,-5,-1,-1);
    score_t jumpScore(-3);
    GlobalJumpAligner<score_t> aligner(scores,jumpScore);
    JumpAlignmentResult<score_t> result;
    aligner.align(
        seq.begin(),seq.end(),
        ref1.begin(),ref1.end(),
        ref2.begin(),ref2.end(),
        result);

    return result;
}


BOOST_AUTO_TEST_CASE( test_GlobalJumpAligner0 )
{
    static const std::string seq("ABABACDCDC");
    static const std::string ref1("ABABA");
    static const std::string ref2("CDCDC");

    JumpAlignmentResult<score_t> result = testAlign(seq,ref1,ref2);

    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align1.apath),"5M");
    BOOST_REQUIRE_EQUAL(result.align1.alignStart,0u);
    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align2.apath),"5M");
    BOOST_REQUIRE_EQUAL(result.align2.alignStart,0u);
}


BOOST_AUTO_TEST_CASE( test_GlobalJumpAligner1 )
{
    static const std::string seq("ABABACDCDC");
    static const std::string ref1("ABABAX");
    static const std::string ref2("CDCDC");

    JumpAlignmentResult<score_t> result = testAlign(seq,ref1,ref2);

    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align1.apath),"5M");
    BOOST_REQUIRE_EQUAL(result.align1.alignStart,0u);
    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align2.apath),"5M");
    BOOST_REQUIRE_EQUAL(result.align2.alignStart,0u);
}


BOOST_AUTO_TEST_CASE( test_GlobalJumpAligner2 )
{
    static const std::string seq("ABABACDCDC");
    static const std::string ref1("ABABA");
    static const std::string ref2("XCDCDC");

    JumpAlignmentResult<score_t> result = testAlign(seq,ref1,ref2);

    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align1.apath),"5M");
    BOOST_REQUIRE_EQUAL(result.align1.alignStart,0u);
    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align2.apath),"5M");
    BOOST_REQUIRE_EQUAL(result.align2.alignStart,1u);
}


BOOST_AUTO_TEST_CASE( test_GlobalJumpAlignerLong )
{
    static const std::string seq("ABABACDCDC");
    static const std::string ref1("dslfjfkjaslABABAlsjfkdsflsk");
    static const std::string ref2("sdfldsklkjdCDCDCfsdlkjfslk");

    JumpAlignmentResult<score_t> result = testAlign(seq,ref1,ref2);

    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align1.apath),"5M");
    BOOST_REQUIRE_EQUAL(result.align1.alignStart,11u);
    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align2.apath),"5M");
    BOOST_REQUIRE_EQUAL(result.align2.alignStart,11u);
}


BOOST_AUTO_TEST_CASE( test_GlobalJumpAlignerSimpleIndels )
{
    static const std::string seq("ABABAABABACDCDCDyCDCDC");
    static const std::string ref1("xABABABABABAx");
    static const std::string ref2("xCDCDCDCDCDCDCx");

    JumpAlignmentResult<score_t> result = testAlign(seq,ref1,ref2);

    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align1.apath),"5M1D5M");
    BOOST_REQUIRE_EQUAL(result.align1.alignStart,1u);
    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align2.apath),"6M1I5M");
    BOOST_REQUIRE_EQUAL(result.align2.alignStart,1u);
}



BOOST_AUTO_TEST_CASE( test_GlobalJumpAlignerBPInsert )
{
    static const std::string seq("ABABABABABA1234CDCDCDCDCDC");
    static const std::string ref1("xABABABABABAx");
    static const std::string ref2("xCDCDCDCDCDCDCx");

    JumpAlignmentResult<score_t> result = testAlign(seq,ref1,ref2);

    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align1.apath),"11M");
    BOOST_REQUIRE_EQUAL(result.align1.alignStart,1u);
    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align2.apath),"11M");
    BOOST_REQUIRE_EQUAL(result.align2.alignStart,1u);
    BOOST_REQUIRE_EQUAL(result.jumpInsertSize,4u);
}


// define behavior when the breakpoint solution is a 1d range
BOOST_AUTO_TEST_CASE( test_GlobalJumpAlignerBPRange )
{
    static const std::string seq("xyzxyzxyzABCABCABCxyzxyzxyz");
    static const std::string ref1("xyzxyzxyzxyzABCABCstustu");
    static const std::string ref2("stustuABCABCxyzxyzxyzxyz");

    JumpAlignmentResult<score_t> result = testAlign(seq,ref1,ref2);

    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align1.apath),"12M");
    BOOST_REQUIRE_EQUAL(result.align1.alignStart,3u);
    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align2.apath),"15M");
    BOOST_REQUIRE_EQUAL(result.align2.alignStart,6u);
}


// define behavior when the breakpoint solution is a 1d range
BOOST_AUTO_TEST_CASE( test_GlobalJumpAlignerBPRange2 )
{
    static const std::string seq("xyzxyzxyzABCABCABCxyzxyzxyz");
    static const std::string ref1("xyzxyzxyzxyzABCABCABCABCABC");
    static const std::string ref2("ABCABCABCABCABCxyzxyzxyzxyz");

    JumpAlignmentResult<score_t> result = testAlign(seq,ref1,ref2);

    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align1.apath),"9M");
    BOOST_REQUIRE_EQUAL(result.align1.alignStart,3u);
    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align2.apath),"18M");
    BOOST_REQUIRE_EQUAL(result.align2.alignStart,6u);
}


BOOST_AUTO_TEST_CASE( test_GlobalJumpAlignerOnly1 )
{
    static const std::string seq("ABABA");
    static const std::string ref1("xABABAx");
    static const std::string ref2("xCDCDCx");

    JumpAlignmentResult<score_t> result = testAlign(seq,ref1,ref2);

    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align1.apath),"5M");
    BOOST_REQUIRE_EQUAL(result.align1.alignStart,1u);
    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align2.apath),"");
    BOOST_REQUIRE_EQUAL(result.align2.alignStart,0u);
}


BOOST_AUTO_TEST_CASE( test_GlobalJumpAlignerOnly2 )
{
    static const std::string seq("CDCDC");
    static const std::string ref1("xABABAx");
    static const std::string ref2("xCDCDCx");

    JumpAlignmentResult<score_t> result = testAlign(seq,ref1,ref2);

    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align1.apath),"");
    BOOST_REQUIRE_EQUAL(result.align1.alignStart,0u);
    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align2.apath),"5M");
    BOOST_REQUIRE_EQUAL(result.align2.alignStart,1u);
}



BOOST_AUTO_TEST_CASE( test_GlobalJumpAlignerOffEdge )
{
    static const std::string seq("123456ABABACDCDC123456");
    static const std::string ref1("xABABAx");
    static const std::string ref2("xCDCDCx");

    JumpAlignmentResult<score_t> result = testAlign(seq,ref1,ref2);

    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align1.apath),"5S6M");
    BOOST_REQUIRE_EQUAL(result.align1.alignStart,0u);
    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align2.apath),"6M5S");
    BOOST_REQUIRE_EQUAL(result.align2.alignStart,1u);
}


BOOST_AUTO_TEST_SUITE_END()

