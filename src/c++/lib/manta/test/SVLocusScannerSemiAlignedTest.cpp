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

#include "htsapi/SimpleAlignment_bam_util.hpp"
#include "htsapi/align_path_bam_util.hpp"
#include "manta/SVLocusScannerSemiAligned.cpp"

#include "boost/scoped_array.hpp"

BOOST_AUTO_TEST_SUITE(test_SVLocusScannerSemiAligned)

/// the reference sequence is assumed to be a clip starting at the query alignment start position:
///
static void semiAlignTestCase(
    const pos_t alignPos,
    const char* querySeq,
    const char* refSeq,
    unsigned&   leadingLength,
    pos_t&      leadingRefPos,
    unsigned&   trailingLength,
    pos_t&      trailingRefPos)
{
  leadingLength  = 0;
  leadingRefPos  = 0;
  trailingLength = 0;
  trailingRefPos = 0;

  const unsigned querySize(strlen(querySeq));

  ALIGNPATH::path_t inputPath;
  inputPath.push_back(ALIGNPATH::path_segment(ALIGNPATH::MATCH, querySize));

  bam_record bamRead;
  bam1_t*    bamDataPtr(bamRead.get_data());
  edit_bam_cigar(inputPath, *bamDataPtr);

  reference_contig_segment testRefSeg;
  testRefSeg.seq() = refSeq;
  testRefSeg.set_offset(alignPos);

  // initialize test qual array to all Q30's:
  boost::scoped_array<uint8_t> qual(new uint8_t[querySize]);
  for (unsigned i(0); i < querySize; ++i) {
    qual[i] = 30;
  }

  edit_bam_read_and_quality(querySeq, qual.get(), *bamDataPtr);

  SimpleAlignment align(getAlignment(bamRead));
  align.pos = alignPos;

  edgePoorAlignmentLength(
      align,
      bamRead.get_bam_read(),
      testRefSeg,
      5,
      leadingLength,
      leadingRefPos,
      trailingLength,
      trailingRefPos);
}

BOOST_AUTO_TEST_CASE(test_getSVCandidatesFromSemiAligned_null)
{
  BOOST_TEST_MESSAGE("SDS MANTA-668");

  static const pos_t alignPos(500);
  static const char  querySeq[] = "AACCCACAAACATCACACACAAGAGTCCAGAGACCGACTTTTTTCTAAAA";
  static const char  refSeq[]   = "AACCCACAAACATCACACACAAGAGTCCAGAGACCGACTTTTTTCTAAAA";

  unsigned leadingLength(0), trailingLength(0);
  pos_t    leadingRefPos(0), trailingRefPos(0);

  semiAlignTestCase(alignPos, querySeq, refSeq, leadingLength, leadingRefPos, trailingLength, trailingRefPos);

  BOOST_REQUIRE_EQUAL(leadingLength, 0u);
  BOOST_REQUIRE_EQUAL(trailingLength, 0u);
  BOOST_REQUIRE_EQUAL(leadingRefPos, alignPos);
  BOOST_REQUIRE_EQUAL(trailingRefPos, alignPos + 50);
}

BOOST_AUTO_TEST_CASE(test_getSVCandidatesFromSemiAligned_leading)
{
  BOOST_TEST_MESSAGE("SDS MANTA-668");

  static const pos_t alignPos(500);
  static const char  querySeq[] = "AACCCACAAACATCACACACAAGAGTCCAGAGACCGACTTTTTTCTAAAA";
  static const char  refSeq[]   = "AACCTTTTTTCATCACACACAAGAGTCCAGAGACCGACTTTTTTCTAAAA";

  unsigned leadingLength(0), trailingLength(0);
  pos_t    leadingRefPos(0), trailingRefPos(0);

  semiAlignTestCase(alignPos, querySeq, refSeq, leadingLength, leadingRefPos, trailingLength, trailingRefPos);

  BOOST_REQUIRE_EQUAL(leadingLength, 10u);
  BOOST_REQUIRE_EQUAL(trailingLength, 0u);
  BOOST_REQUIRE_EQUAL(leadingRefPos, alignPos + 10);
  BOOST_REQUIRE_EQUAL(trailingRefPos, alignPos + 50);
}

BOOST_AUTO_TEST_CASE(test_getSVCandidatesFromSemiAligned_trailing)
{
  BOOST_TEST_MESSAGE("SDS MANTA-668");

  static const pos_t alignPos(500);
  static const char  querySeq[] = "AACCCACAAACATCACACACAAGAGTCCAGAGACCGACTTTTTTCTAAAA";
  static const char  refSeq[]   = "AACCCACAAACATCACACACAAGAGTCCAGAGACCGACTTCCCCCCAAAA";

  unsigned leadingLength(0), trailingLength(0);
  pos_t    leadingRefPos(0), trailingRefPos(0);

  semiAlignTestCase(alignPos, querySeq, refSeq, leadingLength, leadingRefPos, trailingLength, trailingRefPos);

  BOOST_REQUIRE_EQUAL(leadingLength, 0u);
  BOOST_REQUIRE_EQUAL(trailingLength, 10u);
  BOOST_REQUIRE_EQUAL(leadingRefPos, alignPos);
  BOOST_REQUIRE_EQUAL(trailingRefPos, alignPos + 50 - 10);
}

BOOST_AUTO_TEST_CASE(test_getSVCandidatesFromSemiAligned_both)
{
  BOOST_TEST_MESSAGE("SDS MANTA-668");

  static const pos_t alignPos(500);
  static const char  querySeq[] = "AACCCACAAACATCACACACAAGAGTCCAGAGACCGACTTTTTTCTAAAA";
  static const char  refSeq[]   = "AACCTTTTTTCATCACACACAAGAGTCCAGAGACCGACTTCCCCCCAAAA";

  unsigned leadingLength(0), trailingLength(0);
  pos_t    leadingRefPos(0), trailingRefPos(0);

  semiAlignTestCase(alignPos, querySeq, refSeq, leadingLength, leadingRefPos, trailingLength, trailingRefPos);

  BOOST_REQUIRE_EQUAL(leadingLength, 10u);
  BOOST_REQUIRE_EQUAL(trailingLength, 10u);
  BOOST_REQUIRE_EQUAL(leadingRefPos, alignPos + 10);
  BOOST_REQUIRE_EQUAL(trailingRefPos, alignPos + 50 - 10);
}

BOOST_AUTO_TEST_CASE(test_getSVCandidatesFromSemiAligned_mismatch)
{
  BOOST_TEST_MESSAGE("SDS MANTA-668");

  static const pos_t alignPos(500);
  static const char  querySeq[] = "AACCCACAAACATCACACACAAGAGTCCAGAGACCGACTTTTTTCTAAAA";
  static const char  refSeq[]   = "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG";

  unsigned leadingLength(0), trailingLength(0);
  pos_t    leadingRefPos(0), trailingRefPos(0);

  semiAlignTestCase(alignPos, querySeq, refSeq, leadingLength, leadingRefPos, trailingLength, trailingRefPos);

  BOOST_REQUIRE_EQUAL(leadingLength, 50u);
  BOOST_REQUIRE_EQUAL(trailingLength, 50u);
  BOOST_REQUIRE_EQUAL(leadingRefPos, alignPos + 50);
  BOOST_REQUIRE_EQUAL(trailingRefPos, (alignPos + 50 - 50));
}

BOOST_AUTO_TEST_SUITE_END()
