// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Manta - Structural Variant and Indel Caller
// Copyright (c) 2013-2016 Illumina, Inc.
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



BOOST_AUTO_TEST_CASE( test_supportingReadConsistency )
{
    // test against observed case where a single bad read could kill the whole assembly

    SmallAssemblerOptions assembleOpt;

    assembleOpt.minWordLength = 6;
    assembleOpt.maxWordLength = 6;
    assembleOpt.minCoverage = 2;
    assembleOpt.minSeedReads = 3;

    AssemblyReadInput reads;
    reads.push_back(        "AAACGTGTATTA");
    reads.push_back(          "ACGTGTATTACC");
    reads.push_back(           "CGTGTATTACCT");
    reads.push_back(            "GTGTATTACCTA");
    reads.push_back(                "ATTACCTAGTAC");
    reads.push_back(                  "TACCTAGTACTC");
    // the above reads build a contig ACGTG TATTACC TAGTAC
    //
    // Notice ACGTG should not be extended by adding 'A' to the left => AACGTG
    // using the reads below, because they have a different suffix after ACGTG *GCC*
    // Instead, the reads below build a contig CTTA GCTA ACGTG GCC
    reads.push_back("CCCTTAGCTAAC");
    reads.push_back(  "CTTAGCTAACGT");
    reads.push_back(    "TAGCTAACGTGG");
    reads.push_back(      "GCTAACGTGGCC");
    reads.push_back(         "AACGTGGCCTAG");


    AssemblyReadOutput readInfo;
    Assembly contigs;

    runSmallAssembler(assembleOpt, reads, readInfo, contigs);

    BOOST_REQUIRE_EQUAL(contigs.size(),2u);
    BOOST_REQUIRE_EQUAL(contigs[0].seq,"AACGTGTATTACCTAGTAC");
    BOOST_REQUIRE_EQUAL(contigs[1].seq,"CTTAGCTAACGTGGCC");
    for (unsigned i(0); i<6; ++i)
    {
        BOOST_REQUIRE(readInfo[i].isUsed);
        BOOST_REQUIRE_EQUAL(readInfo[i].contigId,0u);
    }

    for (unsigned i(6); i<11; ++i)
    {
        BOOST_REQUIRE(readInfo[i].isUsed);
        BOOST_REQUIRE_EQUAL(readInfo[i].contigId,1u);
    }
}

BOOST_AUTO_TEST_SUITE_END()

