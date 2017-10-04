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

#include "testUtil.hh"
#include "htsapi/SimpleAlignment_bam_util.hh"
#include "manta/SVLocusScanner.cpp"

BOOST_AUTO_TEST_SUITE( test_SVLocusScanner )

// Test all of the paths possible for core read filtering.
BOOST_AUTO_TEST_CASE( test_isReadFilteredCore ) {

    BOOST_TEST_MESSAGE("SDS MANTA-663");

    ALIGNPATH::path_t inputPath;
    cigar_to_apath("100M2000D100M",inputPath);

    bam_record bamRead;
    bam1_t* bamDataPtr(bamRead.get_data());
    edit_bam_cigar(inputPath,*bamDataPtr);

    bam_record primaryFilterRead(bamRead);
    primaryFilterRead.toggle_is_filtered();

    bam_record unmappedRead(bamRead);
    unmappedRead.toggle_is_unmapped();

    bam_record duplicateRead(bamRead);
    duplicateRead.toggle_is_duplicate();

    bam_record supplementNoSASplitRead(bamRead);
    supplementNoSASplitRead.toggle_is_supplementary();
    BOOST_REQUIRE(! supplementNoSASplitRead.isSASplit());

    bam_record supplementSASplitRead(bamRead);
    supplementSASplitRead.toggle_is_supplementary();
    addSupplementaryEvidence(supplementSASplitRead);
    BOOST_REQUIRE(supplementSASplitRead.isSASplit());

    bam_record secondarySASplitRead(bamRead);
    secondarySASplitRead.toggle_is_secondary();
    addSupplementaryEvidence(secondarySASplitRead);
    BOOST_REQUIRE(secondarySASplitRead.isSASplit());

    bam_record secondaryNoSASplitRead(bamRead);
    secondaryNoSASplitRead.toggle_is_secondary();
    BOOST_REQUIRE(! secondaryNoSASplitRead.isSASplit());

    // Test multiple read states for isReadFilteredCore.
    BOOST_REQUIRE(! SVLocusScanner::isReadFilteredCore(bamRead));
    BOOST_REQUIRE(! SVLocusScanner::isReadFilteredCore(unmappedRead));

    BOOST_REQUIRE(SVLocusScanner::isReadFilteredCore(duplicateRead));
    BOOST_REQUIRE(SVLocusScanner::isReadFilteredCore(primaryFilterRead));
    BOOST_REQUIRE(! SVLocusScanner::isReadFilteredCore(secondarySASplitRead));
    BOOST_REQUIRE(SVLocusScanner::isReadFilteredCore(secondaryNoSASplitRead));
    BOOST_REQUIRE(! SVLocusScanner::isReadFilteredCore(supplementSASplitRead));
    BOOST_REQUIRE(SVLocusScanner::isReadFilteredCore(supplementNoSASplitRead));
}

// Test the mapped core read filtering.
BOOST_AUTO_TEST_CASE( test_isMappedReadFilteredCore ) {
    BOOST_TEST_MESSAGE("SDS MANTA-663");

    bam_record bamRead;
    buildBamRecord(bamRead);
    BOOST_REQUIRE(! SVLocusScanner::isMappedReadFilteredCore(bamRead));

    bamRead.toggle_is_unmapped();
    BOOST_REQUIRE(SVLocusScanner::isMappedReadFilteredCore(bamRead));
}

// Test isProperPair for filtering out useless reads.
BOOST_AUTO_TEST_CASE( test_isAnomalousReadPair  ) {
    BOOST_TEST_MESSAGE("SDS MANTA-665");

    bam_header_info bamHeader = bam_header_info();
    bamHeader.chrom_to_index.insert(std::pair<std::string,int32_t>("chrM",0));
    bamHeader.chrom_data.emplace_back("chrM",1000000);
    std::unique_ptr<SVLocusScanner> scanner(buildSVLocusScanner(bamHeader));

    // Proper pair read
    bam_record properPairRead;
    buildBamRecord(properPairRead);
    BOOST_REQUIRE(! scanner->isAnomalousReadPair(properPairRead, 0));

    // Not paired
    bam_record notPairRead;
    buildBamRecord(notPairRead);
    notPairRead.toggle_is_paired();
    BOOST_REQUIRE(scanner->isAnomalousReadPair(notPairRead, 0));

    // Mate not mapped
    bam_record mateNotMappedRead;
    buildBamRecord(mateNotMappedRead);
    mateNotMappedRead .toggle_is_mate_unmapped();
    BOOST_REQUIRE(scanner->isAnomalousReadPair(mateNotMappedRead, 0));

    // small size tests
    // template length = fragment size 0.01
    bam_record lengthSmallRead1;
    buildBamRecord(lengthSmallRead1, 0, 100, 0, 200);
    changeTemplateSize(lengthSmallRead1, 50);
    BOOST_REQUIRE(! scanner->isAnomalousReadPair(lengthSmallRead1, 0));

    // template length < fragment size 0.01 by 1
    bam_record lengthSmallRead2;
    buildBamRecord(lengthSmallRead2, 0, 100, 0, 200);
    changeTemplateSize(lengthSmallRead2, 49);
    BOOST_REQUIRE(scanner->isAnomalousReadPair(lengthSmallRead2, 0));

    // large size tests
    // template length = fragment size 0.99
    bam_record lengthLargeRead1;
    buildBamRecord(lengthLargeRead1, 0, 100, 0, 200);
    changeTemplateSize(lengthLargeRead1, 125);
    BOOST_REQUIRE(! scanner->isAnomalousReadPair(lengthLargeRead1, 0));

    // Proper pair and template length < 0.99 * 1.5 by 1 bp
    bam_record lengthLargeRead2;
    buildBamRecord(lengthLargeRead2, 0, 100, 0, 200);
    changeTemplateSize(lengthLargeRead2, 187);
    BOOST_REQUIRE(! scanner->isAnomalousReadPair(lengthLargeRead2, 0));

    // template length > 0.99 * 1.5
    bam_record lengthLargeRead3;
    buildBamRecord(lengthLargeRead3, 0, 100, 0, 200);
    changeTemplateSize(lengthLargeRead3, 188);
    BOOST_REQUIRE(scanner->isAnomalousReadPair(lengthLargeRead3, 0));
}

// Test the isProperPair functionality for innie-pairs.
BOOST_AUTO_TEST_CASE( test_InniePair  ) {
    BOOST_TEST_MESSAGE("SDS MANTA-665");

    bam_header_info bamHeader = bam_header_info();
    bamHeader.chrom_to_index.insert(std::pair<std::string,int32_t>("chrM",0));
    bamHeader.chrom_data.emplace_back("chrM",1000000);
    std::unique_ptr<SVLocusScanner> scanner(buildSVLocusScanner(bamHeader));

    // Passing condition
    bam_record inniePair;
    buildBamRecord(inniePair);
    BOOST_REQUIRE(! scanner->isAnomalousReadPair(inniePair, 0));

    // different refID
    bam_record difRefIdPair;
    buildBamRecord(difRefIdPair, 0, 100, 1, 200);
    BOOST_REQUIRE(scanner->isAnomalousReadPair(difRefIdPair, 0));

    // same strands (set mate to fwd strand as well.
    bam_record sameStrandInniePair;
    buildBamRecord(sameStrandInniePair);
    sameStrandInniePair.toggle_is_mate_fwd_strand();
    BOOST_REQUIRE(scanner->isAnomalousReadPair(sameStrandInniePair, 0));

    // fwd strand pos > mate pos
    bam_record badPosInniePair;
    buildBamRecord(badPosInniePair, 0, 500, 0, 100);
    sameStrandInniePair.toggle_is_mate_fwd_strand();
    BOOST_REQUIRE(scanner->isAnomalousReadPair(badPosInniePair, 0));

    // pos = mate pos with all strand configurations
    bam_record samePosInniePair;
    buildBamRecord(samePosInniePair, 0, 100, 0, 100);
    BOOST_REQUIRE(! scanner->isAnomalousReadPair(samePosInniePair, 0));
    samePosInniePair.toggle_is_fwd_strand();
    BOOST_REQUIRE(scanner->isAnomalousReadPair(samePosInniePair, 0));
    samePosInniePair.toggle_is_mate_fwd_strand();
    BOOST_REQUIRE(! scanner->isAnomalousReadPair(samePosInniePair, 0));
    samePosInniePair.toggle_is_mate_fwd_strand();
    BOOST_REQUIRE(scanner->isAnomalousReadPair(samePosInniePair, 0));
}

// Test NonCompressedAnomalous read filtering.s
BOOST_AUTO_TEST_CASE( test_LargePair  )
{
    BOOST_TEST_MESSAGE("SDS MANTA-665");

    bam_header_info bamHeader = bam_header_info();
    bamHeader.chrom_to_index.insert(std::pair<std::string,int32_t>("chrM",0));
    bamHeader.chrom_data.emplace_back("chrM",1000000);
    std::unique_ptr<SVLocusScanner> scanner(buildSVLocusScanner(bamHeader));

    // Control - test normal, proper pair against isNonCompressedAnomalousReadPair(..)
    bam_record properPair;
    buildBamRecord(properPair);
    BOOST_REQUIRE(! scanner->isNonCompressedAnomalousReadPair(properPair, 0));

    // refId != mate refID
    bam_record difRefIdPair;
    buildBamRecord(difRefIdPair, 0, 100, 1, 200);
    BOOST_REQUIRE(scanner->isNonCompressedAnomalousReadPair(difRefIdPair, 0));

    // template length = 0
    bam_record lengthZeroRead;
    buildBamRecord(lengthZeroRead, 0, 100, 0, 100);
    changeTemplateSize(lengthZeroRead, 0);
    BOOST_REQUIRE(scanner->isNonCompressedAnomalousReadPair(lengthZeroRead, 0));

    // template length = fragment size 0.99
    bam_record lengthLargeRead1;
    buildBamRecord(lengthLargeRead1, 0, 100, 0, 200);
    changeTemplateSize(lengthLargeRead1, 125);
    BOOST_REQUIRE(! scanner->isNonCompressedAnomalousReadPair(lengthLargeRead1, 0));

    // Proper pair and template length > fragment size 0.99, but less than fragment size 0.99 * 1.5
    bam_record lengthLargeRead2;
    buildBamRecord(lengthLargeRead2, 0, 100, 0, 200);
    changeTemplateSize(lengthLargeRead2, 126);
    BOOST_REQUIRE(! scanner->isNonCompressedAnomalousReadPair(lengthLargeRead2, 0));

    // Not proper pair and template length > fragment size 0.99, but less than template length 0.99 * 1.5
    bam_record lengthLargeRead5;
    buildBamRecord(lengthLargeRead5, 0, 200, 0, 100);
    changeTemplateSize(lengthLargeRead5, 126);
    BOOST_REQUIRE(scanner->isNonCompressedAnomalousReadPair(lengthLargeRead5, 0));

    // template length > 0.99 * 1.5
    bam_record lengthLargeRead3;
    buildBamRecord(lengthLargeRead3, 0, 100, 0, 200);
    changeTemplateSize(lengthLargeRead3, 188);
    BOOST_REQUIRE(scanner->isNonCompressedAnomalousReadPair(lengthLargeRead3, 0));

    // Proper pair and template length < 0.99 * 1.5 by 1 bp
    bam_record lengthLargeRead4;
    buildBamRecord(lengthLargeRead4, 0, 100, 0, 200);
    changeTemplateSize(lengthLargeRead4, 187);
    BOOST_REQUIRE(! scanner->isNonCompressedAnomalousReadPair(lengthLargeRead4, 0));
}

// Test isSVEvidence for filtering out useless reads.
BOOST_AUTO_TEST_CASE( test_isSVEvidence )
{
    BOOST_TEST_MESSAGE("SDS MANTA-666");

    bam_header_info bamHeader = bam_header_info();
    bamHeader.chrom_to_index.insert(std::pair<std::string,int32_t>("chrM",0));
    bamHeader.chrom_data.emplace_back("chrM",1000000);
    std::unique_ptr<SVLocusScanner> scanner(buildSVLocusScanner(bamHeader));

    const unsigned defaultReadGroupIndex = 0;
    reference_contig_segment seq = reference_contig_segment();
    const reference_contig_segment& refSeq(seq);
    SVLocusEvidenceCount* incountsPtr = new SVLocusEvidenceCount();

    // Normal read is not evidence
    bam_record normalRead;
    buildBamRecord(normalRead);
    BOOST_REQUIRE(! scanner->isSVEvidence(normalRead, defaultReadGroupIndex, refSeq, incountsPtr));

    bam_record supplementSASplitRead;
    buildBamRecord(supplementSASplitRead);
    addSupplementaryEvidence(supplementSASplitRead);

    // same read with SA returns true for SV evidence.
    BOOST_REQUIRE(scanner->isSVEvidence(supplementSASplitRead, defaultReadGroupIndex, refSeq, incountsPtr));
}

// Test that different size indels are handled properly.
BOOST_AUTO_TEST_CASE( test_LargeIndels ) {
    BOOST_TEST_MESSAGE("SDS MANTA-667");

    // minCandidateVariantSize = 8
    bam_header_info bamHeader = bam_header_info();
    bamHeader.chrom_to_index.insert(std::pair<std::string,int32_t>("chrM",1000000));
    bamHeader.chrom_data.emplace_back("chrM",1000000);
    std::unique_ptr<SVLocusScanner> scanner(buildSVLocusScanner(bamHeader));

    const unsigned defaultReadGroupIndex = 0;
    reference_contig_segment seq = reference_contig_segment();
    const reference_contig_segment& refSeq(seq);
    SVLocusEvidenceCount* incountsPtr = new SVLocusEvidenceCount();

    // Normal read is not evidence
    bam_record normalRead;
    buildBamRecord(normalRead);
    BOOST_REQUIRE(! scanner->isSVEvidence(normalRead, defaultReadGroupIndex, refSeq, incountsPtr));

    // Non-indel read is not evidence
    bam_record largeSoftclipRead;
    buildBamRecord(largeSoftclipRead, 0, 100, 0, 200, 100, 15, "100M2000S");
    BOOST_REQUIRE(! scanner->isSVEvidence(largeSoftclipRead, defaultReadGroupIndex, refSeq, incountsPtr));

    // large deletion is SV evidence
    bam_record largeDeletionRead;
    buildBamRecord(largeDeletionRead, 0, 100, 0, 200, 100, 15, "100M2000D100M");
    BOOST_REQUIRE(scanner->isSVEvidence(largeDeletionRead, defaultReadGroupIndex, refSeq, incountsPtr));

    // small deletion is not SV evidence
    bam_record smallDeletionRead;
    buildBamRecord(smallDeletionRead, 0, 100, 0, 200, 100, 15, "100M7D100M");
    BOOST_REQUIRE(! scanner->isSVEvidence(smallDeletionRead, defaultReadGroupIndex, refSeq, incountsPtr));

    // Deletion size = minCandidateVariantSize - deletion is SV evidence
    bam_record equalDeletionRead;
    buildBamRecord(equalDeletionRead, 0, 100, 0, 200, 100, 15, "100M8D100M");
    BOOST_REQUIRE(scanner->isSVEvidence(equalDeletionRead, defaultReadGroupIndex, refSeq, incountsPtr));

    // large insertion is SV evidence
    bam_record largeInsertionRead;
    buildBamRecord(largeInsertionRead, 0, 100, 0, 200, 100, 15, "100M2000I100M");
    BOOST_REQUIRE(scanner->isSVEvidence(largeInsertionRead, defaultReadGroupIndex, refSeq, incountsPtr));

    // small insertion is not SV evidence
    bam_record smallInsertionRead;
    buildBamRecord(smallInsertionRead, 0, 100, 0, 200, 100, 15, "100M7I100M");
    BOOST_REQUIRE(! scanner->isSVEvidence(smallInsertionRead, defaultReadGroupIndex, refSeq, incountsPtr));

    // insertion size = minCandidateVariantSize - insertion is SV evidence
    bam_record equalInsertionRead;
    buildBamRecord(equalInsertionRead, 0, 100, 0, 200, 100, 15, "100M8I100M");
    BOOST_REQUIRE(scanner->isSVEvidence(equalInsertionRead, defaultReadGroupIndex, refSeq, incountsPtr));
}

// Test filtering for semiAligned reads.
BOOST_AUTO_TEST_CASE( test_SemiAlignedReads ) {
    BOOST_TEST_MESSAGE("SDS MANTA-668");

    static const pos_t alignPos(500);

    bam_header_info bamHeader = bam_header_info();
    bamHeader.chrom_to_index.insert(std::pair<std::string,int32_t>("chrM",0));
    bamHeader.chrom_data.emplace_back("chrM",1000000);
    std::unique_ptr<SVLocusScanner> scanner(buildSVLocusScanner(bamHeader));

    reference_contig_segment ref = reference_contig_segment();
    static const char refSeq[] =       "AACCTTTTTTCATCACACACAAGAGTCCAGAGACCGACTTCCCCCCAAAA";
    ref.seq() = refSeq;
    ref.set_offset(alignPos);

    // Semi-aligned read - is evidence
    static const char semiQuerySeq[] = "AACCCACAAACATCACACACAAGAGTCCAGAGACCGACTTTTTTCTAAAA";
    static const char smallSeq[] =     "AACCCACAA";
    bam_record semiRead;
    buildBamRecord(semiRead, 0, alignPos, 0, alignPos+50, 50, 15, "50M", semiQuerySeq);
    semiRead.toggle_is_paired();

    SimpleAlignment semiAlignment(getAlignment(semiRead));
    BOOST_REQUIRE(scanner->isSemiAlignedEvidence(semiRead, semiAlignment, ref));

    // normal read - not semi-aligned evidence
    bam_record adaptorPairRead;
    buildBamRecord(adaptorPairRead, 0, alignPos, 0, alignPos+50, 50, 15, "50M", refSeq);
    SimpleAlignment adaptorPairAlignment(getAlignment(adaptorPairRead));
    BOOST_REQUIRE(! scanner->isSemiAlignedEvidence(adaptorPairRead, adaptorPairAlignment, ref));

    // Tests regarding minCandidateVariantSize <= 100
    std::unique_ptr<SVLocusScanner> minCandidate100scanner(buildSVLocusScanner(bamHeader, "tempStatsFile.txt", "tempAlignFile.txt", false, 100));
    BOOST_REQUIRE(minCandidate100scanner->isSVEvidence(semiRead, 0, ref));

    std::unique_ptr<SVLocusScanner> minCandidate101scanner(buildSVLocusScanner(bamHeader, "tempStatsFile.txt", "tempAlignFile.txt", false, 101));
    BOOST_REQUIRE(! minCandidate101scanner->isSVEvidence(semiRead, 0, ref));

    // Test the semi aligned read test. Size and spacing to test is_possible_adapter_pair
    bam_record unmappedRead;
    buildBamRecord(unmappedRead, 0, alignPos, 0, alignPos+5, 50, 15, "50M", semiQuerySeq);
    unmappedRead.toggle_is_unmapped();
    BOOST_REQUIRE(scanner->isSemiAlignedEvidence(unmappedRead, semiAlignment, ref));

    bam_record mateUnmappedRead;
    buildBamRecord(mateUnmappedRead, 0, alignPos, 0, alignPos+50, 50, 15, "50M", semiQuerySeq);
    mateUnmappedRead.toggle_is_mate_unmapped();
    BOOST_REQUIRE(scanner->isSemiAlignedEvidence(mateUnmappedRead, semiAlignment, ref));

    // semiRead different chromosomes
    bam_record diffChromRead;
    buildBamRecord(diffChromRead, 0, alignPos, 1, alignPos+50, 50, 15, "50M", semiQuerySeq);
    BOOST_REQUIRE(scanner->isSemiAlignedEvidence(diffChromRead, semiAlignment, ref));

    // same strand
    bam_record sameStrandRead;
    buildBamRecord(sameStrandRead, 0, alignPos, 1, alignPos+50, 50, 15, "50M", semiQuerySeq);
    sameStrandRead.toggle_is_mate_fwd_strand();
    BOOST_REQUIRE(scanner->isSemiAlignedEvidence(sameStrandRead, semiAlignment, ref));

    // No overlap unless in RNA mode
    bam_record overlappingRead;
    buildBamRecord(overlappingRead, 0, alignPos, 0, alignPos+20, 50, 15, "50M", semiQuerySeq);
    BOOST_REQUIRE(! scanner->isSemiAlignedEvidence(overlappingRead, semiAlignment, ref));

    // Make scanner in RNA mode.
    std::unique_ptr<SVLocusScanner> rnaScanner(buildSVLocusScanner(bamHeader, "tempStatsFile.txt", "tempAlignFile.txt", true, 100));
    BOOST_REQUIRE(rnaScanner->isSemiAlignedEvidence(overlappingRead, semiAlignment, ref));

    // unmapped semi read positions 10bp apart - true
    bam_record pos10bpApartRead;
    buildBamRecord(pos10bpApartRead, 0, alignPos, 0, alignPos+10, 5, 15, "5M", semiQuerySeq);
    SimpleAlignment smallAlignment(getAlignment(pos10bpApartRead));

    //pos10bpApartRead.toggle_is_mate_unmapped();s
    BOOST_REQUIRE(rnaScanner->isSemiAlignedEvidence(pos10bpApartRead, semiAlignment, ref));

    // unmapped semi read positions 9bp apart - false
    bam_record pos9bpApartRead;
    buildBamRecord(pos9bpApartRead, 0, alignPos, 0, alignPos+9, 9, 15, "9M", smallSeq);
    BOOST_REQUIRE(! rnaScanner->isSemiAlignedEvidence(pos9bpApartRead, semiAlignment, ref));

    // unmapped semi read positions -50bp apart - true
    bam_record pos50bpApartRead;
    buildBamRecord(pos50bpApartRead, 0, alignPos, 0, alignPos-50, 20, 15, "20M", semiQuerySeq);
    //pos50bpApartRead.toggle_is_mate_unmapped();
    BOOST_REQUIRE(rnaScanner->isSemiAlignedEvidence(pos50bpApartRead, semiAlignment, ref));

    // unmapped semi read positions -49bp apart - false
    bam_record pos49bpApartRead;
    buildBamRecord(pos49bpApartRead, 0, alignPos, 0, alignPos-49, 20, 15, "20M", semiQuerySeq);
    //mateUnmappedRead.toggle_is_mate_unmapped();
    BOOST_REQUIRE(! rnaScanner->isSemiAlignedEvidence(pos49bpApartRead, semiAlignment, ref));
}

BOOST_AUTO_TEST_CASE( test_getSVCandidatesFromReadIndels )
{
    const bool isRNA(false);
    const bool isStranded(true);
    const ReadScannerOptions opt;
    const ReadScannerDerivOptions dopt(opt,isRNA,isStranded);

    ALIGNPATH::path_t inputPath;
    cigar_to_apath("100M2000D100M",inputPath);

    bam_record bamRead;
    bam1_t* bamDataPtr(bamRead.get_data());
    edit_bam_cigar(inputPath,*bamDataPtr);

    SimpleAlignment align(getAlignment(bamRead));
    align.tid = 0;

    std::vector<SVObservation> candidates;

    bam_header_info hdr_info;
    hdr_info.chrom_data.emplace_back("chr1", 1000000);

    getSVCandidatesFromReadIndels(opt, dopt, align, hdr_info, SourceOfSVEvidenceInDNAFragment::UNKNOWN, candidates);

    BOOST_REQUIRE_EQUAL(candidates.size(),1u);
    BOOST_REQUIRE(candidates[0].bp1.interval.range.is_pos_intersect(100));
    BOOST_REQUIRE(candidates[0].bp2.interval.range.is_pos_intersect(2100));
}

BOOST_AUTO_TEST_SUITE_END()
