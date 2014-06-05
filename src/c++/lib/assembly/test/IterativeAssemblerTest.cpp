// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Manta
// Copyright (c) 2013-2014 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/sequencing/licenses/>
//

///
/// \author Xiaoyu Chen
///

#include "boost/test/unit_test.hpp"

#include "IterativeAssembler.cpp"


BOOST_AUTO_TEST_SUITE( test_IterativeAssembler )

/*
BOOST_AUTO_TEST_CASE( test_CircleDetector )
{
	IterativeAssemblerOptions assembleOpt;
	str_uint_map_t wordCount;
	std::set<std::string> repeatWords;

	wordCount["TACCA"] = 3;
	wordCount["CCACC"] = 3;
	wordCount["CACCA"] = 3;
	wordCount["ACCAC"] = 3;
	wordCount["CCACA"] = 3;
	wordCount["CACAC"] = 3;
	wordCount["ACACA"] = 3;
	wordCount["AAAAA"] = 2;

	getRepeatKmers(assembleOpt, wordCount, repeatWords);

	// the first circle
	BOOST_REQUIRE_EQUAL(repeatWords.count("ACCAC"), 1u);
	BOOST_REQUIRE_EQUAL(repeatWords.count("CACCA"), 1u);
	BOOST_REQUIRE_EQUAL(repeatWords.count("CCACC"), 1u);

	BOOST_REQUIRE_EQUAL(repeatWords.count("TACCA"), 0u);
	BOOST_REQUIRE_EQUAL(repeatWords.count("CCACA"), 0u);

	// the second circle
	BOOST_REQUIRE_EQUAL(repeatWords.count("CACAC"), 1u);
	BOOST_REQUIRE_EQUAL(repeatWords.count("ACACA"), 1u);

	// homopolymer: self-circle
	BOOST_REQUIRE_EQUAL(repeatWords.count("AAAAA"), 1u);
}


BOOST_AUTO_TEST_CASE( test_BasicAssembler )
{
    // test simple assembly functions at a single word size:
    IterativeAssemblerOptions assembleOpt;

    assembleOpt.minWordLength = 6;
    assembleOpt.maxWordLength = 6;

    AssemblyReadInput reads;

    reads.push_back("ACGTGTATTACC");
    reads.push_back(  "GTGTATTACCTA");
    reads.push_back(      "ATTACCTAGTAC");
    reads.push_back(        "TACCTAGTACTC");
    reads.push_back("123456789123");

    AssemblyReadOutput readInfo;
    Assembly contigs;

    runIterativeAssembler(assembleOpt, reads, readInfo, contigs);

    BOOST_REQUIRE_EQUAL(contigs.size(),1u);
    BOOST_REQUIRE_EQUAL(contigs[0].seq,"GTGTATTACCTAGTAC");
    for (unsigned i(0); i<4; ++i)
    {
        BOOST_REQUIRE(readInfo[i].isUsed);
        BOOST_REQUIRE_EQUAL(readInfo[i].contigId,0u);
    }
    BOOST_REQUIRE(! readInfo[4].isUsed);
}*/


BOOST_AUTO_TEST_CASE( test_IterativeKmer )
{
	// test simple assembly functions at a single word size:
	IterativeAssemblerOptions assembleOpt;

	assembleOpt.minWordLength = 3;
	//assembleOpt.maxWordLength = 3;
	assembleOpt.maxWordLength = 7;
	assembleOpt.wordStepSize = 4;
	assembleOpt.minCoverage = 1;
	assembleOpt.minSupportReads = 1;

	AssemblyReadInput reads;

	reads.push_back("ACACACGCCT");
	reads.push_back(      "GCCTTCTCTC");
	reads.push_back("123456789123");

	AssemblyReadOutput readInfo;
	Assembly contigs;

	runIterativeAssembler(assembleOpt, reads, readInfo, contigs);

	BOOST_REQUIRE_EQUAL(contigs.size(),1u);
	BOOST_REQUIRE_EQUAL(contigs[0].seq,"ACACACGCCTTCTCTC");
	for (unsigned i(0); i<2; ++i)
	{
		BOOST_REQUIRE(readInfo[i].isUsed);
		BOOST_REQUIRE_EQUAL(readInfo[i].contigId,0u);
	}
	BOOST_REQUIRE(! readInfo[3].isUsed);
}


BOOST_AUTO_TEST_SUITE_END()

