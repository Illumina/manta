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
// <https://github.com/sequencing/licenses/>
//

#include "boost/test/unit_test.hpp"

#include "blt_util/align_path.hh"


BOOST_AUTO_TEST_SUITE( test_align_path )


BOOST_AUTO_TEST_CASE( test_apath_clean_seqmatch )
{
    const std::string testCigar("10M1D10=2X10=1D1M1=1=1X1=1X");
    ALIGNPATH::path_t testPath;
    cigar_to_apath(testCigar.c_str(), testPath);

    apath_clean_seqmatch(testPath);

    BOOST_REQUIRE_EQUAL(apath_to_cigar(testPath), "10M1D22M1D6M");
}


BOOST_AUTO_TEST_CASE( test_apath_add_seqmatch )
{
    static const std::string testRead("AABAXXXY");
    static const std::string testRef ("AAAADXXXX");

    static const std::string testCigar("4M1D4M");
    ALIGNPATH::path_t testPath;
    cigar_to_apath(testCigar.c_str(), testPath);

    apath_add_seqmatch(
        testRead.begin(), testRead.end(),
        testRef.begin(), testRef.end(),
        testPath);

    BOOST_REQUIRE_EQUAL(apath_to_cigar(testPath), "2=1X1=1D3=1X");
}


BOOST_AUTO_TEST_SUITE_END()

