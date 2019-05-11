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

#include "htsapi/SimpleAlignment_bam_util.hpp"
#include "manta/SVLocusScanner.cpp"
#include "test/testAlignmentDataUtil.hpp"
#include "test/testSVLocusScanner.hpp"
#include "test/testUtil.hpp"

BOOST_AUTO_TEST_SUITE(test_SVLocusScanner)

// Test isProperPair for filtering out useless reads.
BOOST_AUTO_TEST_CASE(test_isAnomalousReadPair)
{
  BOOST_TEST_MESSAGE("SDS MANTA-665");

  const bam_header_info           bamHeader(buildTestBamHeader());
  std::unique_ptr<SVLocusScanner> scanner(buildTestSVLocusScanner(bamHeader));
  std::string                     querySeq = "TCTATCACCCATTTTACCACTCACGGGAGCTCTCC";

  // Proper pair read
  bam_record properPairRead;
  buildTestBamRecord(properPairRead, 0, 100, 0, 200, 35, 20, "35M", querySeq, 100);
  BOOST_REQUIRE(!scanner->isAnomalousReadPair(properPairRead, 0));

  // Not paired
  bam_record notPairRead;
  buildTestBamRecord(notPairRead);
  notPairRead.toggle_is_paired();
  BOOST_REQUIRE(scanner->isAnomalousReadPair(notPairRead, 0));

  // Mate not mapped
  bam_record mateNotMappedRead;
  buildTestBamRecord(mateNotMappedRead);
  mateNotMappedRead.toggle_is_mate_unmapped();
  BOOST_REQUIRE(scanner->isAnomalousReadPair(mateNotMappedRead, 0));

  // small size tests
  // template length = fragment size 0.01
  bam_record lengthSmallRead1;
  buildTestBamRecord(lengthSmallRead1, 0, 100, 0, 200);
  changeTemplateSize(lengthSmallRead1, 50);
  BOOST_REQUIRE(!scanner->isAnomalousReadPair(lengthSmallRead1, 0));

  // template length < fragment size 0.01 by 1
  bam_record lengthSmallRead2;
  buildTestBamRecord(lengthSmallRead2, 0, 100, 0, 200);
  changeTemplateSize(lengthSmallRead2, 49);
  BOOST_REQUIRE(scanner->isAnomalousReadPair(lengthSmallRead2, 0));

  // large size tests
  // template length = fragment size 0.99
  bam_record lengthLargeRead1;
  buildTestBamRecord(lengthLargeRead1, 0, 100, 0, 200);
  changeTemplateSize(lengthLargeRead1, 125);
  BOOST_REQUIRE(!scanner->isAnomalousReadPair(lengthLargeRead1, 0));

  // Proper pair and template length < 0.99 * 1.5 by 1 bp
  bam_record lengthLargeRead2;
  buildTestBamRecord(lengthLargeRead2, 0, 100, 0, 200);
  changeTemplateSize(lengthLargeRead2, 187);
  BOOST_REQUIRE(!scanner->isAnomalousReadPair(lengthLargeRead2, 0));

  // template length > 0.99 * 1.5
  bam_record lengthLargeRead3;
  buildTestBamRecord(lengthLargeRead3, 0, 100, 0, 200);
  changeTemplateSize(lengthLargeRead3, 188);
  BOOST_REQUIRE(scanner->isAnomalousReadPair(lengthLargeRead3, 0));
}

// Test the isProperPair functionality for innie-pairs.
BOOST_AUTO_TEST_CASE(test_InniePair)
{
  BOOST_TEST_MESSAGE("SDS MANTA-665");

  const bam_header_info           bamHeader(buildTestBamHeader());
  std::unique_ptr<SVLocusScanner> scanner(buildTestSVLocusScanner(bamHeader));
  std::string                     querySeq = "TCTATCACCCATTTTACCACTCACGGGAGCTCTCC";

  // Passing condition
  bam_record inniePair;
  buildTestBamRecord(inniePair, 0, 100, 0, 200, 35, 20, "35M", querySeq, 100);
  BOOST_REQUIRE(!scanner->isAnomalousReadPair(inniePair, 0));

  // different refID
  bam_record difRefIdPair;
  buildTestBamRecord(difRefIdPair, 0, 100, 1, 200);
  BOOST_REQUIRE(scanner->isAnomalousReadPair(difRefIdPair, 0));

  // same strands (set mate to fwd strand as well.
  bam_record sameStrandInniePair;
  buildTestBamRecord(sameStrandInniePair);
  sameStrandInniePair.toggle_is_mate_fwd_strand();
  BOOST_REQUIRE(scanner->isAnomalousReadPair(sameStrandInniePair, 0));

  // fwd strand pos > mate pos
  bam_record badPosInniePair;
  buildTestBamRecord(badPosInniePair, 0, 500, 0, 100);
  sameStrandInniePair.toggle_is_mate_fwd_strand();
  BOOST_REQUIRE(scanner->isAnomalousReadPair(badPosInniePair, 0));

  // pos = mate pos with all strand configurations
  bam_record samePosInniePair;
  buildTestBamRecord(samePosInniePair, 0, 100, 0, 100, 35, 20, "35M", querySeq, 100);
  BOOST_REQUIRE(!scanner->isAnomalousReadPair(samePosInniePair, 0));
  samePosInniePair.toggle_is_fwd_strand();
  BOOST_REQUIRE(scanner->isAnomalousReadPair(samePosInniePair, 0));
  samePosInniePair.toggle_is_mate_fwd_strand();
  BOOST_REQUIRE(!scanner->isAnomalousReadPair(samePosInniePair, 0));
  samePosInniePair.toggle_is_mate_fwd_strand();
  BOOST_REQUIRE(scanner->isAnomalousReadPair(samePosInniePair, 0));
}

// Test NonCompressedAnomalous read filtering.s
BOOST_AUTO_TEST_CASE(test_LargePair)
{
  BOOST_TEST_MESSAGE("SDS MANTA-665");

  const bam_header_info           bamHeader(buildTestBamHeader());
  std::unique_ptr<SVLocusScanner> scanner(buildTestSVLocusScanner(bamHeader));
  std::string                     querySeq = "TCTATCACCCATTTTACCACTCACGGGAGCTCTCC";

  // Control - test normal, proper pair against isNonCompressedAnomalousReadPair(..)
  bam_record properPair;
  buildTestBamRecord(properPair, 0, 100, 0, 200, 35, 20, "35M", querySeq, 100);
  BOOST_REQUIRE(!scanner->isNonCompressedAnomalousReadPair(properPair, 0));

  // refId != mate refID
  bam_record difRefIdPair;
  buildTestBamRecord(difRefIdPair, 0, 100, 1, 200);
  BOOST_REQUIRE(scanner->isNonCompressedAnomalousReadPair(difRefIdPair, 0));

  // template length = 0
  bam_record lengthZeroRead;
  buildTestBamRecord(lengthZeroRead, 0, 100, 0, 100);
  changeTemplateSize(lengthZeroRead, 0);
  BOOST_REQUIRE(scanner->isNonCompressedAnomalousReadPair(lengthZeroRead, 0));

  // template length = fragment size 0.99
  bam_record lengthLargeRead1;
  buildTestBamRecord(lengthLargeRead1, 0, 100, 0, 200);
  changeTemplateSize(lengthLargeRead1, 125);
  BOOST_REQUIRE(!scanner->isNonCompressedAnomalousReadPair(lengthLargeRead1, 0));

  // Proper pair and template length > fragment size 0.99, but less than fragment size 0.99 * 1.5
  bam_record lengthLargeRead2;
  buildTestBamRecord(lengthLargeRead2, 0, 100, 0, 200);
  changeTemplateSize(lengthLargeRead2, 126);
  BOOST_REQUIRE(!scanner->isNonCompressedAnomalousReadPair(lengthLargeRead2, 0));

  // Not proper pair and template length > fragment size 0.99, but less than template length 0.99 * 1.5
  bam_record lengthLargeRead5;
  buildTestBamRecord(lengthLargeRead5, 0, 200, 0, 100);
  changeTemplateSize(lengthLargeRead5, 126);
  BOOST_REQUIRE(scanner->isNonCompressedAnomalousReadPair(lengthLargeRead5, 0));

  // template length > 0.99 * 1.5
  bam_record lengthLargeRead3;
  buildTestBamRecord(lengthLargeRead3, 0, 100, 0, 200);
  changeTemplateSize(lengthLargeRead3, 188);
  BOOST_REQUIRE(scanner->isNonCompressedAnomalousReadPair(lengthLargeRead3, 0));

  // Proper pair and template length < 0.99 * 1.5 by 1 bp
  bam_record lengthLargeRead4;
  buildTestBamRecord(lengthLargeRead4, 0, 100, 0, 200);
  changeTemplateSize(lengthLargeRead4, 187);
  BOOST_REQUIRE(!scanner->isNonCompressedAnomalousReadPair(lengthLargeRead4, 0));
}

// Test isSVEvidence for filtering out useless reads.
BOOST_AUTO_TEST_CASE(test_isSVEvidence)
{
  BOOST_TEST_MESSAGE("SDS MANTA-666");

  const bam_header_info           bamHeader(buildTestBamHeader());
  std::unique_ptr<SVLocusScanner> scanner(buildTestSVLocusScanner(bamHeader));
  std::string                     querySeq = "TCTATCACCCATTTTACCACTCACGGGAGCTCTCC";

  const unsigned                        defaultReadGroupIndex = 0;
  reference_contig_segment              seq                   = reference_contig_segment();
  const reference_contig_segment&       refSeq(seq);
  std::unique_ptr<SVLocusEvidenceCount> incountsPtr(new SVLocusEvidenceCount());

  // Normal read is not evidence
  bam_record normalRead;
  buildTestBamRecord(normalRead, 0, 100, 0, 200, 35, 20, "35M", querySeq, 100);
  BOOST_REQUIRE(!scanner->isSVEvidence(normalRead, defaultReadGroupIndex, refSeq, incountsPtr.get()));

  bam_record supplementSASplitRead;
  buildTestBamRecord(supplementSASplitRead);
  addSupplementaryAlignmentEvidence(supplementSASplitRead);

  // same read with SA returns true for SV evidence.
  BOOST_REQUIRE(
      scanner->isSVEvidence(supplementSASplitRead, defaultReadGroupIndex, refSeq, incountsPtr.get()));
}

// Test that different size indels are handled properly.
BOOST_AUTO_TEST_CASE(test_LargeIndels)
{
  BOOST_TEST_MESSAGE("SDS MANTA-667");

  // minCandidateVariantSize = 8
  const bam_header_info           bamHeader(buildTestBamHeader());
  std::unique_ptr<SVLocusScanner> scanner(buildTestSVLocusScanner(bamHeader));
  std::string                     querySeq = "TCTATCACCCATTTTACCACTCACGGGAGCTCTCC";

  const unsigned                        defaultReadGroupIndex = 0;
  reference_contig_segment              seq                   = reference_contig_segment();
  const reference_contig_segment&       refSeq(seq);
  std::unique_ptr<SVLocusEvidenceCount> incountsPtr(new SVLocusEvidenceCount());

  // Normal read is not evidence
  bam_record normalRead;
  buildTestBamRecord(normalRead, 0, 100, 0, 200, 35, 20, "35M", querySeq, 100);
  BOOST_REQUIRE(!scanner->isSVEvidence(normalRead, defaultReadGroupIndex, refSeq, incountsPtr.get()));

  // Non-indel read is not evidence
  bam_record largeSoftclipRead;
  buildTestBamRecord(largeSoftclipRead, 0, 100, 0, 200, 2100, 15, "100M2000S", "", 100);
  BOOST_REQUIRE(!scanner->isSVEvidence(largeSoftclipRead, defaultReadGroupIndex, refSeq, incountsPtr.get()));

  // large deletion is SV evidence
  bam_record largeDeletionRead;
  buildTestBamRecord(largeDeletionRead, 0, 100, 0, 200, 200, 15, "100M2000D100M");
  BOOST_REQUIRE(scanner->isSVEvidence(largeDeletionRead, defaultReadGroupIndex, refSeq, incountsPtr.get()));

  // small deletion is not SV evidence
  bam_record smallDeletionRead;
  buildTestBamRecord(smallDeletionRead, 0, 100, 0, 200, 50, 15, "25M7D25M", "", 100);
  BOOST_REQUIRE(!scanner->isSVEvidence(smallDeletionRead, defaultReadGroupIndex, refSeq, incountsPtr.get()));

  // Deletion size = minCandidateVariantSize - deletion is SV evidence
  bam_record equalDeletionRead;
  buildTestBamRecord(equalDeletionRead, 0, 100, 0, 200, 200, 15, "100M8D100M");
  BOOST_REQUIRE(scanner->isSVEvidence(equalDeletionRead, defaultReadGroupIndex, refSeq, incountsPtr.get()));

  // large insertion is SV evidence
  bam_record largeInsertionRead;
  buildTestBamRecord(largeInsertionRead, 0, 100, 0, 200, 2200, 15, "100M2000I100M");
  BOOST_REQUIRE(scanner->isSVEvidence(largeInsertionRead, defaultReadGroupIndex, refSeq, incountsPtr.get()));

  // small insertion is not SV evidence
  bam_record smallInsertionRead;
  buildTestBamRecord(smallInsertionRead, 0, 100, 0, 200, 27, 15, "10M7I10M", "", 100);
  BOOST_REQUIRE(!scanner->isSVEvidence(smallInsertionRead, defaultReadGroupIndex, refSeq, incountsPtr.get()));

  // insertion size = minCandidateVariantSize - insertion is SV evidence
  bam_record equalInsertionRead;
  buildTestBamRecord(equalInsertionRead, 0, 100, 0, 200, 208, 15, "100M8I100M");
  BOOST_REQUIRE(scanner->isSVEvidence(equalInsertionRead, defaultReadGroupIndex, refSeq, incountsPtr.get()));
}

// Test filtering for semiAligned reads.
BOOST_AUTO_TEST_CASE(test_SemiAlignedReads)
{
  BOOST_TEST_MESSAGE("SDS MANTA-668");

  static const pos_t alignPos(500);

  const bam_header_info           bamHeader(buildTestBamHeader());
  std::unique_ptr<SVLocusScanner> scanner(buildTestSVLocusScanner(bamHeader));

  reference_contig_segment ref      = reference_contig_segment();
  static const char        refSeq[] = "AACCTTTTTTCATCACACACAAGAGTCCAGAGACCGACTTCCCCCCAAAA";
  ref.seq()                         = refSeq;
  ref.set_offset(alignPos);

  // Semi-aligned read - is evidence
  static const char semiQuerySeq[] = "AACCCACAAACATCACACACAAGAGTCCAGAGACCGACTTTTTTCTAAAA";
  bam_record        semiRead;
  buildTestBamRecord(semiRead, 0, alignPos, 0, alignPos + 50, 50, 15, "50M", semiQuerySeq);
  semiRead.toggle_is_paired();

  SimpleAlignment semiAlignment(getAlignment(semiRead));
  BOOST_REQUIRE(scanner->isSemiAlignedEvidence(semiRead, semiAlignment, ref));

  // normal read - not semi-aligned evidence
  bam_record adaptorPairRead;
  buildTestBamRecord(adaptorPairRead, 0, alignPos, 0, alignPos + 50, 50, 15, "50M", refSeq);
  SimpleAlignment adaptorPairAlignment(getAlignment(adaptorPairRead));
  BOOST_REQUIRE(!scanner->isSemiAlignedEvidence(adaptorPairRead, adaptorPairAlignment, ref));

  // Tests regarding minCandidateVariantSize <= 100
  std::unique_ptr<SVLocusScanner> minCandidate100scanner(buildTestSVLocusScanner(bamHeader, false, 100));
  BOOST_REQUIRE(minCandidate100scanner->isSVEvidence(semiRead, 0, ref));

  std::unique_ptr<SVLocusScanner> minCandidate101scanner(buildTestSVLocusScanner(bamHeader, false, 101));
  BOOST_REQUIRE(!minCandidate101scanner->isSVEvidence(semiRead, 0, ref));

  // Test the semi aligned read test. Size and spacing to test is_possible_adapter_pair
  bam_record unmappedRead;
  buildTestBamRecord(unmappedRead, 0, alignPos, 0, alignPos + 5, 50, 15, "50M", semiQuerySeq);
  unmappedRead.toggle_is_unmapped();
  BOOST_REQUIRE(scanner->isSemiAlignedEvidence(unmappedRead, semiAlignment, ref));

  bam_record mateUnmappedRead;
  buildTestBamRecord(mateUnmappedRead, 0, alignPos, 0, alignPos + 50, 50, 15, "50M", semiQuerySeq);
  mateUnmappedRead.toggle_is_mate_unmapped();
  BOOST_REQUIRE(scanner->isSemiAlignedEvidence(mateUnmappedRead, semiAlignment, ref));

  // semiRead different chromosomes
  bam_record diffChromRead;
  buildTestBamRecord(diffChromRead, 0, alignPos, 1, alignPos + 50, 50, 15, "50M", semiQuerySeq);
  BOOST_REQUIRE(scanner->isSemiAlignedEvidence(diffChromRead, semiAlignment, ref));

  // same strand
  bam_record sameStrandRead;
  buildTestBamRecord(sameStrandRead, 0, alignPos, 1, alignPos + 50, 50, 15, "50M", semiQuerySeq);
  sameStrandRead.toggle_is_mate_fwd_strand();
  BOOST_REQUIRE(scanner->isSemiAlignedEvidence(sameStrandRead, semiAlignment, ref));

  // No overlap unless in RNA mode
  bam_record overlappingRead;
  buildTestBamRecord(overlappingRead, 0, alignPos, 0, alignPos + 20, 50, 15, "50M", semiQuerySeq);
  BOOST_REQUIRE(!scanner->isSemiAlignedEvidence(overlappingRead, semiAlignment, ref));

  // Make scanner in RNA mode.
  std::unique_ptr<SVLocusScanner> rnaScanner(buildTestSVLocusScanner(bamHeader, true, 100));
  BOOST_REQUIRE(rnaScanner->isSemiAlignedEvidence(overlappingRead, semiAlignment, ref));
}

BOOST_AUTO_TEST_CASE(test_getSVCandidatesFromReadIndels)
{
  const bool                    isStranded(true);
  const ReadScannerOptions      opt;
  const ReadScannerDerivOptions dopt(opt, isStranded);

  ALIGNPATH::path_t inputPath;
  cigar_to_apath("100M2000D100M", inputPath);

  bam_record bamRead;
  bam1_t*    bamDataPtr(bamRead.get_data());
  edit_bam_cigar(inputPath, *bamDataPtr);

  SimpleAlignment align(getAlignment(bamRead));
  align.tid = 0;

  std::vector<SVObservation> candidates;

  bam_header_info hdr_info;
  hdr_info.chrom_data.emplace_back("chr1", 1000000);

  getSVCandidatesFromReadIndels(
      opt, dopt, align, hdr_info, SourceOfSVEvidenceInDNAFragment::UNKNOWN, candidates);

  BOOST_REQUIRE_EQUAL(candidates.size(), 1u);
  BOOST_REQUIRE(candidates[0].bp1.interval.range.is_pos_intersect(100));
  BOOST_REQUIRE(candidates[0].bp2.interval.range.is_pos_intersect(2100));
}

BOOST_AUTO_TEST_SUITE_END()
