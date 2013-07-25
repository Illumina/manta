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
    AlignmentScores<score_t> scores(2,-4,-5,-1);
    score_t jumpScore(-4);
    GlobalJumpAligner<score_t> aligner(scores,jumpScore);
    JumpAlignmentResult<score_t> result;
    aligner.align(
            seq.begin(),seq.end(),
            ref1.begin(),ref1.end(),
            ref2.begin(),ref2.end(),
            result);

    return result;
}


BOOST_AUTO_TEST_CASE( test_GlobalJumpALigner1 )
{
    static const std::string seq("ABABACDCDC");
    static const std::string ref1("ABABABABAB");
    static const std::string ref2("XCDCDCDCDCD");

    JumpAlignmentResult<score_t> result = testAlign(seq,ref1,ref2);

    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align2.apath),"5M");
    BOOST_REQUIRE_EQUAL(result.align2.alignStart,1);
    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align1.apath),"5M");
    BOOST_REQUIRE_EQUAL(result.align1.alignStart,0);
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

