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
/// \author Chris Saunders
///

#include "boost/test/unit_test.hpp"

#include "SmallAssembler.hh"


BOOST_AUTO_TEST_SUITE( test_SmallAssembler )


BOOST_AUTO_TEST_CASE( test_SmallAssembler1 )
{
    // test simple assembly functions at a single word size:

    SmallAssemblerOptions assembleOpt;

    assembleOpt.minWordLength = 6;
    assembleOpt.maxWordLength = 6;
    assembleOpt.minCoverage = 2;
    assembleOpt.minSeedReads = 3;

    AssemblyReadInput reads;

    reads.push_back("ACGTGTATTACC");
    reads.push_back(  "GTGTATTACCTA");
    reads.push_back(      "ATTACCTAGTAC");
    reads.push_back(        "TACCTAGTACTC");
    reads.push_back("123456789123");

    AssemblyReadOutput readInfo;
    Assembly contigs;

    runSmallAssembler(assembleOpt, reads, readInfo, contigs);

    BOOST_REQUIRE_EQUAL(contigs.size(),1u);
    BOOST_REQUIRE_EQUAL(contigs[0].seq,"GTGTATTACCTAGTAC");
    for (unsigned i(0); i<4; ++i)
    {
        BOOST_REQUIRE(readInfo[i].isUsed);
        BOOST_REQUIRE_EQUAL(readInfo[i].contigId,0u);
    }
    BOOST_REQUIRE(! readInfo[4].isUsed);
}



BOOST_AUTO_TEST_CASE( test_PoisonRead )
{
    // test against observed case where a single bad read could kill the whole assembly

    SmallAssemblerOptions assembleOpt;

    assembleOpt.minWordLength = 6;
    assembleOpt.maxWordLength = 6;
    assembleOpt.minCoverage = 2;
    assembleOpt.minSeedReads = 3;

    AssemblyReadInput reads;

    reads.push_back("ACGTGTATTACC");
    reads.push_back(  "GTGTATTACCTA");
    reads.push_back(      "ATTACCTAGTAC");
    reads.push_back(        "TACCTAGTACTC");
    reads.push_back("AAAAAAAAAAAAAAAAAAAA");

    AssemblyReadOutput readInfo;
    Assembly contigs;

    runSmallAssembler(assembleOpt, reads, readInfo, contigs);

    BOOST_REQUIRE_EQUAL(contigs.size(),1u);
    BOOST_REQUIRE_EQUAL(contigs[0].seq,"GTGTATTACCTAGTAC");
    for (unsigned i(0); i<4; ++i)
    {
        BOOST_REQUIRE(readInfo[i].isUsed);
        BOOST_REQUIRE_EQUAL(readInfo[i].contigId,0u);
    }
    BOOST_REQUIRE(readInfo[4].isUsed);
    BOOST_REQUIRE_EQUAL(readInfo[4].contigId,0u);
}


BOOST_AUTO_TEST_SUITE_END()

