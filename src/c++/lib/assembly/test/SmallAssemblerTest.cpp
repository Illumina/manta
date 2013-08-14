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
    reads.push_back("GTGTATTACCTA");
    reads.push_back("ATTACCTAGTAC");
    reads.push_back("TACCTAGTACTC");
    reads.push_back("123456789123");

    AssemblyReadOutput readInfo;
    Assembly contigs;

    runSmallAssembler(assembleOpt, reads, readInfo, contigs);
    
    BOOST_REQUIRE_EQUAL(contigs.size(),1u);
    BOOST_REQUIRE_EQUAL(contigs[0].seq,"GTGTATTACCTAGTAC");
    for(unsigned i(0);i<4;++i)
    {
        BOOST_REQUIRE(readInfo[i].isUsed);
        BOOST_REQUIRE_EQUAL(readInfo[i].contigId,0u);
    }
    BOOST_REQUIRE(! readInfo[5].isUsed);


}


BOOST_AUTO_TEST_SUITE_END()

