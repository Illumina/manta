//
// Manta - Structural Variant and Indel Caller
// Copyright (c) 2013-2025 Illumina, Inc.
//
// This program is licensed under the terms of the Polyform strict license
//
// ***As far as the law allows, the software comes as is, without
// any warranty or condition, and the licensor will not be liable
// to you for any damages arising out of these terms or the use
// or nature of the software, under any kind of legal claim.***
//
// You should have received a copy of the PolyForm Strict License 1.0.0
// along with this program.  If not, see <https://polyformproject.org/licenses/strict/1.0.0>.
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
