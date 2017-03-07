// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Manta - Structural Variant and Indel Caller
// Copyright (c) 2013-2017 Illumina, Inc.
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
/// \author Chris Saunders
///

#include "boost/test/unit_test.hpp"

// test static function in TU:
#include "applications/GenerateSVCandidates/SVCandidateAssemblyRefiner.cpp"


BOOST_AUTO_TEST_SUITE( test_SVRefiner )


BOOST_AUTO_TEST_CASE( test_GetVariantRange )
{
    known_pos_range2 res;

    const std::string seq1("ABCDDABC");
    const std::string seq2("ABCDDDABC");

    {
        // left shifted case:
        const known_pos_range2 seq1Range(3,3);
        const known_pos_range2 seq2Range(3,4);

        // order reflects a deletion
        res = getVariantRange(seq2,seq2Range,seq1,seq1Range);
        BOOST_REQUIRE_EQUAL(res.begin_pos(), 0);
        BOOST_REQUIRE_EQUAL(res.end_pos(), 2);

        // order reflects an insertion
        res = getVariantRange(seq1,seq1Range,seq2,seq2Range);
        BOOST_REQUIRE_EQUAL(res.begin_pos(), 0);
        BOOST_REQUIRE_EQUAL(res.end_pos(), 2);
    }

    {
        // right shifted case:
        const known_pos_range2 seq1Range(5,5);
        const known_pos_range2 seq2Range(5,6);

        // order reflects a deletion
        res = getVariantRange(seq2,seq2Range,seq1,seq1Range);
        BOOST_REQUIRE_EQUAL(res.begin_pos(), -2);
        BOOST_REQUIRE_EQUAL(res.end_pos(), 0);

        // order reflects an insertion
        res = getVariantRange(seq1,seq1Range,seq2,seq2Range);
        BOOST_REQUIRE_EQUAL(res.begin_pos(), -2);
        BOOST_REQUIRE_EQUAL(res.end_pos(), 0);
    }
}


BOOST_AUTO_TEST_SUITE_END()
