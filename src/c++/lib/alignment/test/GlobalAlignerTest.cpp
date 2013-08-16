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

#include "GlobalAligner.hh"

#include "blt_util/align_path.hh"

#include <string>



BOOST_AUTO_TEST_SUITE( test_GlobalAligner )

typedef short int score_t;

static
AlignmentResult<score_t>
testAlign(
    const std::string& seq,
    const std::string& ref)
{
    AlignmentScores<score_t> scores(2,-4,-5,-1,-4);
    GlobalAligner<score_t> aligner(scores);
    AlignmentResult<score_t> result;
    aligner.align(seq.begin(),seq.end(),ref.begin(),ref.end(),result);

    return result;
}


BOOST_AUTO_TEST_CASE( test_GlobalAligner1 )
{
    static const std::string seq("D");
    static const std::string ref("ABCDEF");

    AlignmentResult<score_t> result = testAlign(seq,ref);

    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align.apath),"1M");
    BOOST_REQUIRE_EQUAL(result.align.beginPos,3u);
}



BOOST_AUTO_TEST_CASE( test_GlobalAlignerDelete )
{
    static const std::string seq("BCDEFHIKLM");
    static const std::string ref("ABCDEFGHIKLMN");

    AlignmentResult<score_t> result = testAlign(seq,ref);

    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align.apath),"5M1D5M");
    BOOST_REQUIRE_EQUAL(result.align.beginPos,1u);
}


BOOST_AUTO_TEST_CASE( test_GlobalAlignerInsert )
{
    static const std::string seq("BCDEFGXHIKLM");
    static const std::string ref("ABCDEFGHIKLMN");

    AlignmentResult<score_t> result = testAlign(seq,ref);

    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align.apath),"6M1I5M");
    BOOST_REQUIRE_EQUAL(result.align.beginPos,1u);
}


BOOST_AUTO_TEST_CASE( test_GlobalAlignerInsertDelete )
{
    static const std::string seq("BBBBBBCDXYZHIKLMMMM");
    static const std::string ref("ABBBBBBCDEFGHIKLMMMMN");

    AlignmentResult<score_t> result = testAlign(seq,ref);

    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align.apath),"8M3I3D8M");
    BOOST_REQUIRE_EQUAL(result.align.beginPos,1u);
}

BOOST_AUTO_TEST_CASE( test_GlobalAlignerInsertDelete2 )
{
    static const std::string seq("BBBBBBCDEXYHIKLMMMM");
    static const std::string ref("ABBBBBBCDEFGHIKLMMMMN");

    AlignmentResult<score_t> result = testAlign(seq,ref);

    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align.apath),"19M");
    BOOST_REQUIRE_EQUAL(result.align.beginPos,1u);
}


BOOST_AUTO_TEST_CASE( test_GlobalAlignerShortRef1 )
{
    static const std::string seq("ABCD");
    static const std::string ref("BCD");

    AlignmentResult<score_t> result = testAlign(seq,ref);

    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align.apath),"1S3M");
    BOOST_REQUIRE_EQUAL(result.align.beginPos,0u);
    BOOST_REQUIRE_EQUAL(result.score,2);
}


BOOST_AUTO_TEST_CASE( test_GlobalAlignerShortRef2 )
{
    static const std::string seq("ABCD");
    static const std::string ref("ABC");

    AlignmentResult<score_t> result = testAlign(seq,ref);

    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align.apath),"3M1S");
    BOOST_REQUIRE_EQUAL(result.align.beginPos,0u);
    BOOST_REQUIRE_EQUAL(result.score,2);
}


BOOST_AUTO_TEST_CASE( test_GlobalAlignerShortRef3 )
{
    static const std::string seq("ABCD");
    static const std::string ref("B");

    AlignmentResult<score_t> result = testAlign(seq,ref);

    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align.apath),"1S1M2S");
    BOOST_REQUIRE_EQUAL(result.align.beginPos,0u);
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
    BOOST_REQUIRE_EQUAL(result.align.beginPos,0u);
}

// show that the method left aligns an insertion within a repeat
//
BOOST_AUTO_TEST_CASE( test_GlobalAlignerLeftShift2 )
{
    static const std::string seq("ABCDEFFFFFFFGHIJKL");
    static const std::string ref("ABCDEFFFFFFGHIJKL");

    AlignmentResult<score_t> result = testAlign(seq,ref);

    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align.apath),"5M1I12M");
    BOOST_REQUIRE_EQUAL(result.align.beginPos,0u);
}


BOOST_AUTO_TEST_SUITE_END()

