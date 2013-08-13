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
    BOOST_REQUIRE_EQUAL(result.align1.alignStart,0);
    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align2.apath),"5M");
    BOOST_REQUIRE_EQUAL(result.align2.alignStart,0);
}


BOOST_AUTO_TEST_CASE( test_GlobalJumpAligner1 )
{
    static const std::string seq("ABABACDCDC");
    static const std::string ref1("ABABAX");
    static const std::string ref2("CDCDC");

    JumpAlignmentResult<score_t> result = testAlign(seq,ref1,ref2);

    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align1.apath),"5M");
    BOOST_REQUIRE_EQUAL(result.align1.alignStart,0);
    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align2.apath),"5M");
    BOOST_REQUIRE_EQUAL(result.align2.alignStart,0);
}


BOOST_AUTO_TEST_CASE( test_GlobalJumpAligner2 )
{
    static const std::string seq("ABABACDCDC");
    static const std::string ref1("ABABA");
    static const std::string ref2("XCDCDC");

    JumpAlignmentResult<score_t> result = testAlign(seq,ref1,ref2);

    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align1.apath),"5M");
    BOOST_REQUIRE_EQUAL(result.align1.alignStart,0);
    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align2.apath),"5M");
    BOOST_REQUIRE_EQUAL(result.align2.alignStart,1);
}

BOOST_AUTO_TEST_CASE( test_GlobalJumpAlignerLong )
{
    static const std::string seq("ABABACDCDC");
    static const std::string ref1("dslfjfkjaslABABAlsjfkdsflsk");
    static const std::string ref2("sdfldsklkjdCDCDCfsdlkjfslk");

    JumpAlignmentResult<score_t> result = testAlign(seq,ref1,ref2);

    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align1.apath),"5M");
    BOOST_REQUIRE_EQUAL(result.align1.alignStart,11);
    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align2.apath),"5M");
    BOOST_REQUIRE_EQUAL(result.align2.alignStart,11);
}


BOOST_AUTO_TEST_CASE( test_GlobalJumpAlignerSimpleIndels )
{
    static const std::string seq("ABABAABABACDCDCDyCDCDC");
    static const std::string ref1("xABABABABABAx");
    static const std::string ref2("xCDCDCDCDCDCDCx");

    JumpAlignmentResult<score_t> result = testAlign(seq,ref1,ref2);

    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align1.apath),"5M1D5M");
    BOOST_REQUIRE_EQUAL(result.align1.alignStart,1);
    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align2.apath),"6M1I5M");
    BOOST_REQUIRE_EQUAL(result.align2.alignStart,1);
}



BOOST_AUTO_TEST_CASE( test_GlobalJumpAlignerBPInsert )
{
    static const std::string seq("ABABABABABA1234CDCDCDCDCDC");
    static const std::string ref1("xABABABABABAx");
    static const std::string ref2("xCDCDCDCDCDCDCx");

    JumpAlignmentResult<score_t> result = testAlign(seq,ref1,ref2);

    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align1.apath),"11M");
    BOOST_REQUIRE_EQUAL(result.align1.alignStart,1);
    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align2.apath),"11M");
    BOOST_REQUIRE_EQUAL(result.align2.alignStart,1);
    BOOST_REQUIRE_EQUAL(result.jumpInsertSize,4);
}


// define behavior when the breakpoint solution is a 1d range
BOOST_AUTO_TEST_CASE( test_GlobalJumpAlignerBPRange )
{
    static const std::string seq("xyzxyzxyzABCABCABCxyzxyzxyz");
    static const std::string ref1("xyzxyzxyzxyzABCABCstustu");
    static const std::string ref2("stustuABCABCxyzxyzxyzxyz");

    JumpAlignmentResult<score_t> result = testAlign(seq,ref1,ref2);

    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align1.apath),"12M");
    BOOST_REQUIRE_EQUAL(result.align1.alignStart,3);
    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align2.apath),"15M");
    BOOST_REQUIRE_EQUAL(result.align2.alignStart,6);
}


// define behavior when the breakpoint solution is a 1d range
BOOST_AUTO_TEST_CASE( test_GlobalJumpAlignerBPRange2 )
{
    static const std::string seq("xyzxyzxyzABCABCABCxyzxyzxyz");
    static const std::string ref1("xyzxyzxyzxyzABCABCABCABCABC");
    static const std::string ref2("ABCABCABCABCABCxyzxyzxyzxyz");

    JumpAlignmentResult<score_t> result = testAlign(seq,ref1,ref2);

    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align1.apath),"9M");
    BOOST_REQUIRE_EQUAL(result.align1.alignStart,3);
    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align2.apath),"18M");
    BOOST_REQUIRE_EQUAL(result.align2.alignStart,6);
}


BOOST_AUTO_TEST_CASE( test_GlobalJumpAlignerOnly1 )
{
    static const std::string seq("ABABA");
    static const std::string ref1("xABABAx");
    static const std::string ref2("xCDCDCx");

    JumpAlignmentResult<score_t> result = testAlign(seq,ref1,ref2);

    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align1.apath),"5M");
    BOOST_REQUIRE_EQUAL(result.align1.alignStart,1);
    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align2.apath),"");
    BOOST_REQUIRE_EQUAL(result.align2.alignStart,0);
}


BOOST_AUTO_TEST_CASE( test_GlobalJumpAlignerOnly2 )
{
    static const std::string seq("CDCDC");
    static const std::string ref1("xABABAx");
    static const std::string ref2("xCDCDCx");

    JumpAlignmentResult<score_t> result = testAlign(seq,ref1,ref2);

    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align1.apath),"");
    BOOST_REQUIRE_EQUAL(result.align1.alignStart,0);
    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align2.apath),"5M");
    BOOST_REQUIRE_EQUAL(result.align2.alignStart,1);
}



BOOST_AUTO_TEST_CASE( test_GlobalJumpAlignerOffEdge )
{
    static const std::string seq("123456ABABACDCDC123456");
    static const std::string ref1("xABABAx");
    static const std::string ref2("xCDCDCx");

    JumpAlignmentResult<score_t> result = testAlign(seq,ref1,ref2);

    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align1.apath),"5S6M");
    BOOST_REQUIRE_EQUAL(result.align1.alignStart,0);
    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align2.apath),"6M5S");
    BOOST_REQUIRE_EQUAL(result.align2.alignStart,1);
}




#if 0
BOOST_AUTO_TEST_CASE( test_GlobalAlignerDelete )
{
    static const std::string seq("BCDEFHIKLM");
    static const std::string ref("ABCDEFGHIKLMN");

    AlignmentResult<score_t> result = testAlign(seq,ref);

    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align.apath),"5M1D5M");
    BOOST_REQUIRE_EQUAL(result.align.alignStart,1);
}


BOOST_AUTO_TEST_CASE( test_GlobalAlignerInsert )
{
    static const std::string seq("BCDEFGXHIKLM");
    static const std::string ref("ABCDEFGHIKLMN");

    AlignmentResult<score_t> result = testAlign(seq,ref);

    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align.apath),"6M1I5M");
    BOOST_REQUIRE_EQUAL(result.align.alignStart,1);
}


BOOST_AUTO_TEST_CASE( test_GlobalAlignerInsertDelete )
{
    static const std::string seq("BBBBBBCDXYZHIKLMMMM");
    static const std::string ref("ABBBBBBCDEFGHIKLMMMMN");

    AlignmentResult<score_t> result = testAlign(seq,ref);

    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align.apath),"8M3I3D8M");
    BOOST_REQUIRE_EQUAL(result.align.alignStart,1);
}

BOOST_AUTO_TEST_CASE( test_GlobalAlignerInsertDelete2 )
{
    static const std::string seq("BBBBBBCDEXYHIKLMMMM");
    static const std::string ref("ABBBBBBCDEFGHIKLMMMMN");

    AlignmentResult<score_t> result = testAlign(seq,ref);

    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align.apath),"19M");
    BOOST_REQUIRE_EQUAL(result.align.alignStart,1);
}


BOOST_AUTO_TEST_CASE( test_GlobalAlignerShortRef1 )
{
    static const std::string seq("ABCD");
    static const std::string ref("BCD");

    AlignmentResult<score_t> result = testAlign(seq,ref);

    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align.apath),"1S3M");
    BOOST_REQUIRE_EQUAL(result.align.alignStart,0);
    BOOST_REQUIRE_EQUAL(result.score,2);
}


BOOST_AUTO_TEST_CASE( test_GlobalAlignerShortRef2 )
{
    static const std::string seq("ABCD");
    static const std::string ref("ABC");

    AlignmentResult<score_t> result = testAlign(seq,ref);

    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align.apath),"3M1S");
    BOOST_REQUIRE_EQUAL(result.align.alignStart,0);
    BOOST_REQUIRE_EQUAL(result.score,2);
}


BOOST_AUTO_TEST_CASE( test_GlobalAlignerShortRef3 )
{
    static const std::string seq("ABCD");
    static const std::string ref("B");

    AlignmentResult<score_t> result = testAlign(seq,ref);

    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align.apath),"1S1M2S");
    BOOST_REQUIRE_EQUAL(result.align.alignStart,0);
    BOOST_REQUIRE_EQUAL(result.score,-10);
}


// show that the method left aligns a deletion within a repeat
//
BOOST_AUTO_TEST_CASE( test_GlobalAlignerLeftShift )
{
    static const std::string seq("ABCDEFFFFFGHIJKL");
    static const std::string ref("ABCDEFFFFFFGHIJKL");

    AlignmentResult<score_t> result = testAlign(seq,ref);

    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align.apath),"5M1D11M");
    BOOST_REQUIRE_EQUAL(result.align.alignStart,0);
}

// show that the method left aligns an insertion within a repeat
//
BOOST_AUTO_TEST_CASE( test_GlobalAlignerLeftShift2 )
{
    static const std::string seq("ABCDEFFFFFFFGHIJKL");
    static const std::string ref("ABCDEFFFFFFGHIJKL");

    AlignmentResult<score_t> result = testAlign(seq,ref);

    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align.apath),"5M1I12M");
    BOOST_REQUIRE_EQUAL(result.align.alignStart,0);
}
#endif

BOOST_AUTO_TEST_SUITE_END()

