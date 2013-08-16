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

#include "AlignmentUtil.hh"

#include "blt_util/align_path.hh"

#include <string>



BOOST_AUTO_TEST_SUITE( test_AlignmnetUtil )



BOOST_AUTO_TEST_CASE( test_getQuerySegments )
{
    const std::string querySeq("AAAACCGGG");

    JumpAlignmentResult<int> result;

    cigar_to_apath("4M",result.align1.apath);
    result.jumpInsertSize = 2;
    cigar_to_apath("3M",result.align2.apath);

    std::string bp1Seq,insertSeq,bp2Seq;
    getFwdStrandQuerySegments(result,querySeq,false,true,false,
        bp1Seq,bp2Seq,insertSeq);

    BOOST_REQUIRE_EQUAL(bp1Seq,"TTTT");
    BOOST_REQUIRE_EQUAL(insertSeq,"GG");
    BOOST_REQUIRE_EQUAL(bp2Seq,"GGG");
}


BOOST_AUTO_TEST_SUITE_END()

