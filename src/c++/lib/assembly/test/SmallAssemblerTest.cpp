//
// Manta - Structural Variant and Indel Caller
// Copyright (c) 2013-2019 Illumina, Inc.
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

#include "boost/test/unit_test.hpp"

#include "SmallAssembler.hpp"

BOOST_AUTO_TEST_SUITE(test_SmallAssembler)

BOOST_AUTO_TEST_CASE(test_SmallAssembler1)
{
  // test simple assembly functions at a single word size:

  SmallAssemblerOptions assembleOpt;

  assembleOpt.minWordLength = 6;
  assembleOpt.maxWordLength = 6;
  assembleOpt.minCoverage   = 2;
  assembleOpt.minSeedReads  = 3;

  AssemblyReadInput reads;

  // clang-format off
  reads.emplace_back("ACGTGTATTACC");
  reads.emplace_back(  "GTGTATTACCTA");
  reads.emplace_back(      "ATTACCTAGTAC");
  reads.emplace_back(        "TACCTAGTACTC");
  reads.emplace_back("123456789123");
  // clang-format on

  AssemblyReadOutput readInfo;
  Assembly           contigs;

  runSmallAssembler(assembleOpt, reads, readInfo, contigs);

  BOOST_REQUIRE_EQUAL(contigs.size(), 1u);
  BOOST_REQUIRE_EQUAL(contigs[0].seq, "GTGTATTACCTAGTAC");
  for (unsigned i(0); i < 4; ++i) {
    BOOST_REQUIRE(readInfo[i].isUsed);
    BOOST_REQUIRE_EQUAL(readInfo[i].contigIds[0], 0u);
  }
  BOOST_REQUIRE(!readInfo[4].isUsed);
}

BOOST_AUTO_TEST_CASE(test_PoisonRead)
{
  // test against observed case where a single bad read could kill the whole assembly

  SmallAssemblerOptions assembleOpt;

  assembleOpt.minWordLength = 6;
  assembleOpt.maxWordLength = 6;
  assembleOpt.minCoverage   = 2;
  assembleOpt.minSeedReads  = 3;

  AssemblyReadInput reads;

  // clang-format off
  reads.emplace_back("ACGTGTATTACC");
  reads.emplace_back(  "GTGTATTACCTA");
  reads.emplace_back(      "ATTACCTAGTAC");
  reads.emplace_back(        "TACCTAGTACTC");
  reads.emplace_back("AAAAAAAAAAAAAAAAAAAA");
  // clang-format on

  AssemblyReadOutput readInfo;
  Assembly           contigs;

  runSmallAssembler(assembleOpt, reads, readInfo, contigs);

  BOOST_REQUIRE_EQUAL(contigs.size(), 1u);
  BOOST_REQUIRE_EQUAL(contigs[0].seq, "GTGTATTACCTAGTAC");
  for (unsigned i(0); i < 4; ++i) {
    BOOST_REQUIRE(readInfo[i].isUsed);
    BOOST_REQUIRE_EQUAL(readInfo[i].contigIds[0], 0u);
  }
  BOOST_REQUIRE(readInfo[4].isUsed);
  BOOST_REQUIRE_EQUAL(readInfo[4].contigIds.size(), 0u);
}

BOOST_AUTO_TEST_CASE(test_supportingReadConsistency)
{
  // test against observed case where a single bad read could kill the whole assembly

  SmallAssemblerOptions assembleOpt;

  assembleOpt.minWordLength = 6;
  assembleOpt.maxWordLength = 6;
  assembleOpt.minCoverage   = 2;
  assembleOpt.minSeedReads  = 3;

  AssemblyReadInput reads;

  // clang-format off
  reads.emplace_back(        "AAACGTGTATTA");
  reads.emplace_back(          "ACGTGTATTACC");
  reads.emplace_back(           "CGTGTATTACCT");
  reads.emplace_back(            "GTGTATTACCTA");
  reads.emplace_back(                "ATTACCTAGTAC");
  reads.emplace_back(                  "TACCTAGTACTC");
  // clang-format on

  // the above reads build a contig ACGTG TATTACC TAGTAC
  //
  // Notice ACGTG should not be extended by adding 'A' to the left => AACGTG
  // using the reads below, because they have a different suffix after ACGTG *GCC*
  // Instead, the reads below build a contig CTTA GCTA ACGTG GCC

  // clang-format off
  reads.emplace_back("CCCTTAGCTAAC");
  reads.emplace_back(  "CTTAGCTAACGT");
  reads.emplace_back(    "TAGCTAACGTGG");
  reads.emplace_back(      "GCTAACGTGGCC");
  reads.emplace_back(         "AACGTGGCCTAG");
  // clang-format on

  AssemblyReadOutput readInfo;
  Assembly           contigs;

  runSmallAssembler(assembleOpt, reads, readInfo, contigs);

  BOOST_REQUIRE_EQUAL(contigs.size(), 2u);
  BOOST_REQUIRE_EQUAL(contigs[0].seq, "AACGTGTATTACCTAGTAC");
  BOOST_REQUIRE_EQUAL(contigs[1].seq, "CTTAGCTAACGTGGCC");
  for (unsigned i(0); i < 6; ++i) {
    BOOST_REQUIRE(readInfo[i].isUsed);
    BOOST_REQUIRE_EQUAL(readInfo[i].contigIds[0], 0u);
  }

  for (unsigned i(6); i < 11; ++i) {
    BOOST_REQUIRE(readInfo[i].isUsed);
    BOOST_REQUIRE_EQUAL(readInfo[i].contigIds[0], 1u);
  }
}

BOOST_AUTO_TEST_SUITE_END()
