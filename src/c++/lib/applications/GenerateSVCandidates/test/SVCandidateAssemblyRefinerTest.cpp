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
#include "manta/BamStreamerUtils.hpp"
#include "test/testAlignmentDataUtil.hpp"
#include "test/testFileMakers.hpp"
// test static function in TU:
#include "applications/GenerateSVCandidates/SVCandidateAssemblyRefiner.cpp"

BOOST_AUTO_TEST_SUITE(test_SVRefiner)

// Create Temporary bam streams of a bam file which contains
// three interchromosomal read pairs.
struct BamStream {
  BamStream()
  {
    const bam_header_info bamHeader(buildTestBamHeader());

    bam_record bamRecord1;
    buildTestBamRecord(bamRecord1, 0, 8, 1, 69, 35, 15, "35M", "GTCTATCACCCTATTAACCACTCACGGGAGAAAAA", 0);
    bamRecord1.set_qname("Read-1");
    bam_record bamRecord2;
    buildTestBamRecord(bamRecord2, 0, 9, 1, 70, 35, 15, "35M", "TCTATCACCCTATTAACCACTCACGGGAGAAAAAG", 0);
    bamRecord2.set_qname("Read-2");
    bam_record bamRecord3;
    buildTestBamRecord(bamRecord3, 0, 10, 1, 69, 35, 15, "35M", "CTATCACCCTATTAACCACTCACGGGAGAAAAAGA", 0);
    bamRecord3.set_qname("Read-3");
    bam_record bamRecord4;
    buildTestBamRecord(bamRecord4, 0, 299, 0, 350, 35, 15, "35M", "AAACCCCCCCCTCCCCCCGCTTCTGGCCTTTTTTT");
    bamRecord4.set_qname("Read-4");
    bam_record bamRecord5;
    buildTestBamRecord(bamRecord5, 0, 300, 0, 350, 35, 15, "35M", "AACCCCCCCCTCCCCCCGCTTCTGGCCTTTTTTTT");
    bamRecord5.set_qname("Read-5");
    bam_record bamRecord6;
    buildTestBamRecord(bamRecord6, 0, 301, 0, 350, 35, 15, "35M", "ACCCCCCCCTCCCCCCGCTTCTGGCCTTTTTTTTT");
    bamRecord6.set_qname("Read-6");
    bam_record bamRecord7;
    buildTestBamRecord(bamRecord7, 1, 69, 0, 8, 35, 50, "35M", "AAAAAAACTCATCAGTTGATGATACGCCCGAGCAG", 0);
    bamRecord7.toggle_is_mate_fwd_strand();
    bamRecord7.toggle_is_second();
    bamRecord7.toggle_is_fwd_strand();
    bamRecord7.set_qname("Read-1");
    bam_record bamRecord8;
    buildTestBamRecord(bamRecord8, 1, 70, 0, 9, 35, 50, "35M", "AAAAAACTCATCAGTTGATGATACGCCCGAGCAGA", 0);
    bamRecord8.set_qname("Read-2");
    bamRecord8.toggle_is_second();
    bamRecord8.toggle_is_mate_fwd_strand();
    bamRecord8.toggle_is_fwd_strand();
    bam_record bamRecord9;
    buildTestBamRecord(bamRecord9, 1, 71, 0, 10, 35, 50, "35M", "AAAAACTCATCAGTTGATGATACGCCCGAGCAGAT", 0);
    bamRecord9.toggle_is_mate_fwd_strand();
    bamRecord9.toggle_is_second();
    bamRecord9.toggle_is_fwd_strand();
    bamRecord9.set_qname("Read-3");
    readsToAdd.push_back(bamRecord1);
    readsToAdd.push_back(bamRecord2);
    readsToAdd.push_back(bamRecord3);
    readsToAdd.push_back(bamRecord4);
    readsToAdd.push_back(bamRecord5);
    readsToAdd.push_back(bamRecord6);
    readsToAdd.push_back(bamRecord7);
    readsToAdd.push_back(bamRecord8);
    readsToAdd.push_back(bamRecord9);
    bamFileName = _bamFilename();
    buildTestBamFile(bamHeader, readsToAdd, bamFileName);

    const std::string                          referenceFilename = getTestReferenceFilename();
    std::vector<std::string>                   bamFilenames      = {bamFileName};
    std::vector<std::shared_ptr<bam_streamer>> bamStreams;
    openBamStreams(referenceFilename, bamFilenames, bamStreams);
    bamStream = bamStreams[0];
  }

  std::shared_ptr<bam_streamer> bamStream;
  std::vector<bam_record>       readsToAdd;
  std::string                   bamFileName;

private:
  const std::string&     _bamFilename() const { return _bamFilenameMaker.getFilename(); }
  const BamFilenameMaker _bamFilenameMaker;
};

BOOST_FIXTURE_TEST_SUITE(SVRefiner_test_suite, BamStream)

BOOST_AUTO_TEST_CASE(test_GetVariantRange)
{
  known_pos_range2 res;

  const std::string seq1("ABCDDABC");
  const std::string seq2("ABCDDDABC");

  {
    // left shifted case:
    const known_pos_range2 seq1Range(3, 3);
    const known_pos_range2 seq2Range(3, 4);

    // order reflects a deletion
    res = getVariantRange(seq2, seq2Range, seq1, seq1Range);
    BOOST_REQUIRE_EQUAL(res.begin_pos(), 0);
    BOOST_REQUIRE_EQUAL(res.end_pos(), 2);

    // order reflects an insertion
    res = getVariantRange(seq1, seq1Range, seq2, seq2Range);
    BOOST_REQUIRE_EQUAL(res.begin_pos(), 0);
    BOOST_REQUIRE_EQUAL(res.end_pos(), 2);
  }

  {
    // right shifted case:
    const known_pos_range2 seq1Range(5, 5);
    const known_pos_range2 seq2Range(5, 6);

    // order reflects a deletion
    res = getVariantRange(seq2, seq2Range, seq1, seq1Range);
    BOOST_REQUIRE_EQUAL(res.begin_pos(), -2);
    BOOST_REQUIRE_EQUAL(res.end_pos(), 0);

    // order reflects an insertion
    res = getVariantRange(seq1, seq1Range, seq2, seq2Range);
    BOOST_REQUIRE_EQUAL(res.begin_pos(), -2);
    BOOST_REQUIRE_EQUAL(res.end_pos(), 0);
  }
}

// Test start and end coordinate of a read segment based on the cigar elements(like 35M, 5I etc.).
// For example:
// Example-1: Let's say CIGAR of a read is 35M5I15M and segment of cigar is (1, 1), then
//            read segment start and end is 35 and 40 respectively.
// Example-2: Let's say CIGAR of a read is 35M5D5I15M and segment of cigar is (1, 1), then
//            read segment start and end is 35 and 35 respectively.
// Following match types are included in this calculation:
//  1. MATCH
//  2. INSERT
//  3. SOFT_CLIP
//  4. SEQ_MATCH
//  5. SEQ_MISMATCH
BOOST_AUTO_TEST_CASE(test_getInsertTrim)
{
  bam_record bamRecord1;
  buildTestBamRecord(bamRecord1, 0, 9, 0, 90, 55, 15, "35M5I15M");
  SimpleAlignment simpleAlignment1(getAlignment(bamRecord1));
  // segment range is (1, 1)
  known_pos_range2 range1(getInsertTrim(simpleAlignment1.path, std::pair<unsigned, unsigned>(1, 1)));
  // As segment range (1, 1), read segment start will be after 1st element of cigar(35M) and read segment
  // end will be after 2nd element of the cigar. so here read segment start is 35 and read
  // segment end is 40.
  BOOST_REQUIRE_EQUAL(range1.begin_pos(), 35);
  BOOST_REQUIRE_EQUAL(range1.end_pos(), 40);

  bam_record bamRecord2;
  buildTestBamRecord(bamRecord2, 0, 9, 0, 90, 55, 15, "35M5D5I15M");
  SimpleAlignment simpleAlignment2(getAlignment(bamRecord2));
  // segment range is (1, 1)
  known_pos_range2 range2(getInsertTrim(simpleAlignment2.path, std::pair<unsigned, unsigned>(1, 1)));
  // As segment range (1, 1), read segment start will be after 1st element of cigar(35M) and read segment
  // end will be after 2nd element(5D) of the cigar. so here read segment start is 35 and read
  // segment end is 35 (Deletion is not included in the calculation).
  BOOST_REQUIRE_EQUAL(range2.begin_pos(), 35);
  BOOST_REQUIRE_EQUAL(range2.end_pos(), 35);

  bam_record bamRecord3;
  buildTestBamRecord(bamRecord3, 0, 9, 0, 90, 50, 15, "35M5D15M");
  SimpleAlignment simpleAlignment3(getAlignment(bamRecord3));
  // segment range is (1, 2)
  known_pos_range2 range3(getInsertTrim(simpleAlignment3.path, std::pair<unsigned, unsigned>(1, 2)));
  // As segment range (1, 1), read segment start will be after 1st element of cigar(35M) and read
  // end will be after 3rd element(15M) of the cigar. so here read segment start is 35 and read
  // segment end is 50 (Deletion is not included).
  BOOST_REQUIRE_EQUAL(range3.begin_pos(), 35);
  BOOST_REQUIRE_EQUAL(range3.end_pos(), 50);
}

// Test the following cases
// 1. Size of largest indel in a cigar segment
// 2. Size of largest indel with adjacent deletion and insertion
// 3. Size of largest indel among multiple cigar segments
BOOST_AUTO_TEST_CASE(test_getLargestIndelSize)
{
  std::string       testCigar1("35M5I30M");
  ALIGNPATH::path_t path1;
  cigar_to_apath(testCigar1.c_str(), path1);
  std::vector<std::pair<unsigned, unsigned>> segments;
  segments.push_back(std::pair<unsigned, unsigned>(1, 1));
  // Here only segment is (1, 1). As there is only one insertion
  // element in this segment, so size of the largest indel is zero.
  BOOST_REQUIRE_EQUAL(getLargestIndelSize(path1, segments), 5);

  segments.clear();
  segments.push_back(std::pair<unsigned, unsigned>(0, 0));
  segments.push_back(std::pair<unsigned, unsigned>(2, 2));
  // Here segments are (0, 0) and (2, 2). As there is no indel
  // in these two segments of the cigar, so size of the largest
  // indel is zero.
  BOOST_REQUIRE_EQUAL(getLargestIndelSize(path1, segments), 0);

  segments.clear();
  segments.push_back(std::pair<unsigned, unsigned>(0, 0));
  segments.push_back(std::pair<unsigned, unsigned>(1, 1));
  // Here segments are (0, 0) and (1, 1). As there is no indel
  // in the (0, 0) segment and one insertion element(5I) in the (1, 1)
  // segment of the cigar, so size of the largest indel is 5.
  BOOST_REQUIRE_EQUAL(getLargestIndelSize(path1, segments), 5);

  std::string       testCigar2("35M5I30M6D10M");
  ALIGNPATH::path_t path2;
  cigar_to_apath(testCigar2.c_str(), path2);
  segments.clear();
  segments.push_back(std::pair<unsigned, unsigned>(0, 0));
  segments.push_back(std::pair<unsigned, unsigned>(1, 1));
  segments.push_back(std::pair<unsigned, unsigned>(4, 4));
  // Here segments are (0, 0), (1, 1) and (4, 4). Among these 3 segments
  // (1,1) segment has the largest indel and size is 5.
  BOOST_REQUIRE_EQUAL(getLargestIndelSize(path2, segments), 5);

  segments.clear();
  segments.push_back(std::pair<unsigned, unsigned>(0, 0));
  segments.push_back(std::pair<unsigned, unsigned>(1, 1));
  segments.push_back(std::pair<unsigned, unsigned>(3, 4));
  // Here segments are (0, 0), (1, 1) and (3, 4). Among these 3 segments,
  // (3,4) segment has the largest indel and size is 6.
  BOOST_REQUIRE_EQUAL(getLargestIndelSize(path2, segments), 6);

  std::string       testCigar3("35M5I6D40M");
  ALIGNPATH::path_t path3;
  cigar_to_apath(testCigar3.c_str(), path3);
  segments.clear();
  segments.push_back(std::pair<unsigned, unsigned>(1, 2));
  // Here only one segment (1, 2). It has 5I and 6D, so maximum
  // indel size is 6.
  BOOST_REQUIRE_EQUAL(getLargestIndelSize(path3, segments), 6);
  segments.clear();
  segments.push_back(std::pair<unsigned, unsigned>(0, 0));
  segments.push_back(std::pair<unsigned, unsigned>(1, 1));
  segments.push_back(std::pair<unsigned, unsigned>(3, 3));
  // Here segments are (0, 0), (1, 1) and (3, 3). Among these 3 segments
  // (1,1) segment has the largest indel and size is 5.
  BOOST_REQUIRE_EQUAL(getLargestIndelSize(path3, segments), 5);

  segments.clear();
  segments.push_back(std::pair<unsigned, unsigned>(0, 0));
  segments.push_back(std::pair<unsigned, unsigned>(1, 1));
  segments.push_back(std::pair<unsigned, unsigned>(2, 3));
  // Here segments are (0, 0), (1, 1) and (2, 3). Among these 3 segments,
  // (2, 3) segment has the largest indel and size is 6.
  BOOST_REQUIRE_EQUAL(getLargestIndelSize(path3, segments), 6);
}

// 1. Test the larger indel location in the cigar based on
//    minIndel size as indel theshold.
// 2. When there is an adjacent stretch of insertion and deletion,
//    if one of the adjacent segments satisfy minIndel threshold,
//    then whole stretch is treated as a single indel segment.
//    see the example below.
BOOST_AUTO_TEST_CASE(test_getLargeIndelSegments)
{
  std::string       testCigar1("35M5I30M6D10M3I");
  ALIGNPATH::path_t path1;
  cigar_to_apath(testCigar1.c_str(), path1);
  // location of the required indel segments
  std::vector<std::pair<unsigned, unsigned>> segments;

  // There is no indel with size 10 (Indel threshold).
  // So segments size should be zero.
  getLargeIndelSegments(path1, 10, segments);
  BOOST_REQUIRE(segments.size() == 0);

  // Indel threshold = 6. There is only one deletion with size 6.
  // So Segments size should be 1.
  getLargeIndelSegments(path1, 6, segments);
  BOOST_REQUIRE(segments.size() == 1);
  // location of the deletion segment in the cigar.
  BOOST_REQUIRE(segments[0].first == 3);
  BOOST_REQUIRE(segments[0].second == 3);

  // Indel threshold = 3. There are three indels with size greater
  // than or equal to 3. So Segments size should be 3.
  getLargeIndelSegments(path1, 3, segments);
  BOOST_REQUIRE(segments.size() == 3);
  // location of the indel segments.
  BOOST_REQUIRE(segments[0].first == 1);
  BOOST_REQUIRE(segments[0].second == 1);
  BOOST_REQUIRE(segments[1].first == 3);
  BOOST_REQUIRE(segments[1].second == 3);
  BOOST_REQUIRE(segments[2].first == 5);
  BOOST_REQUIRE(segments[2].second == 5);

  segments.clear();
  // Where Insertion and Deletion (5I6D) are adjacent to each other.
  // Indel threshold = 4. Both of the segments (5I & 6D) are satisfied
  // min indel threshold criteria. So Segment start is 1 and segment end is 2.
  std::string       testCigar2("35M5I6D40M3I");
  ALIGNPATH::path_t path2;
  cigar_to_apath(testCigar2.c_str(), path2);
  getLargeIndelSegments(path2, 4, segments);
  BOOST_REQUIRE(segments.size() == 1);
  // location of the indel segments.
  BOOST_REQUIRE(segments[0].first == 1);
  BOOST_REQUIRE(segments[0].second == 2);

  segments.clear();
  // Where Insertion and Deletion (5I3D) are adjacent to each other.
  // Indel threshold = 4. Although Only one segment (5I) is satisfied
  // min indel threshold criteria, still this stretch is treated as single
  // segment (Segment start is 1 and segment end is 2) as one of the segments
  // satisfied the criteria.
  std::string       testCigar3("35M5I3D40M3I");
  ALIGNPATH::path_t path3;
  cigar_to_apath(testCigar3.c_str(), path3);
  getLargeIndelSegments(path3, 4, segments);
  BOOST_REQUIRE(segments.size() == 1);
  // location of the indel segments.
  BOOST_REQUIRE(segments[0].first == 1);
  BOOST_REQUIRE(segments[0].second == 2);

  segments.clear();
  // Where Insertion and Deletion (3I5D) are adjacent to each other.
  // Indel threshold = 4. Although Only one segment (5D) is satisfied
  // min indel threshold criteria, still this stretch is treated as single
  // segment (Segment start is 1 and segment end is 2) as one of the segments
  // satisfied the criteria.
  std::string       testCigar4("35M3I5D40M3I");
  ALIGNPATH::path_t path4;
  cigar_to_apath(testCigar4.c_str(), path4);
  getLargeIndelSegments(path4, 4, segments);
  BOOST_REQUIRE(segments.size() == 1);
  // location of the indel segments.
  BOOST_REQUIRE(segments[0].first == 1);
  BOOST_REQUIRE(segments[0].second == 2);

  segments.clear();
  // Where Insertion and Deletion are adjacent to each other.
  // Indel threshold = 4. None of the segments satisfy indel threshold
  // criteria. So Segments size should be 0.
  std::string       testCigar5("35M3I3D40M3I");
  ALIGNPATH::path_t path5;
  cigar_to_apath(testCigar5.c_str(), path5);
  getLargeIndelSegments(path5, 4, segments);
  BOOST_REQUIRE(segments.size() == 0);
}

// Test the largest insert segment with minInsertThreshold.
BOOST_AUTO_TEST_CASE(test_getLargestInsertSegment)
{
  std::string       testCigar1("35M5I30M6D10M2I");
  ALIGNPATH::path_t path1;
  cigar_to_apath(testCigar1.c_str(), path1);
  // Segments with largest insert segment
  std::vector<std::pair<unsigned, unsigned>> segments1;
  // There is no insertion with size 10 (Insertion threshold).
  // So segments size should be zero.
  getLargestInsertSegment(path1, 10, segments1);
  BOOST_REQUIRE(segments1.size() == 0);
  // Insertion threshold = 3. There is only one insertion segment which
  // is more than 3. So Segments size should be 1.
  getLargestInsertSegment(path1, 3, segments1);
  BOOST_REQUIRE(segments1.size() == 1);
  // Location of the largest insert segment
  BOOST_REQUIRE(segments1[0].first == 1);
  BOOST_REQUIRE(segments1[0].second == 1);

  std::string       testCigar2("35M5I30M6D10M6I");
  ALIGNPATH::path_t path2;
  cigar_to_apath(testCigar2.c_str(), path2);
  std::vector<std::pair<unsigned, unsigned>> segments2;
  // Insertion threshold = 3. There are two insertion segments which
  // are more than 3 but largest insert segment is 6.
  // So Segments size should be 1.
  getLargestInsertSegment(path2, 3, segments2);
  BOOST_REQUIRE(segments2.size() == 1);
  // Location of the largest segment
  BOOST_REQUIRE(segments2[0].first == 5);
  BOOST_REQUIRE(segments2[0].second == 5);

  // Where Insertion and Deletion are adjacent.
  std::string       testCigar3("35M5I6D40M6I");
  ALIGNPATH::path_t path3;
  cigar_to_apath(testCigar3.c_str(), path3);
  std::vector<std::pair<unsigned, unsigned>> segments3;
  // Insertion threshold = 3. There are two insertion segments which
  // are more than 3 but largest insert segment is 6.
  // So Segments size should be 1.
  getLargestInsertSegment(path3, 3, segments3);
  BOOST_REQUIRE(segments3.size() == 1);
  // Location of the largest segment
  BOOST_REQUIRE(segments3[0].first == 4);
  BOOST_REQUIRE(segments3[0].second == 4);
}

// Test the number of locations when query sequence aligns gaplessly with
// the target sequence provided a threshold of mismatch rate(fraction of mismatch
// allowed).
// If target size is 100 and query size is 50, so from 0 to 49th location
// it will try to align and calculate the mismatch rate.
BOOST_AUTO_TEST_CASE(test_getQuerySeqMatchCount)
{
  std::string targetSeq =
      "GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGT"
      "ATTTTCGTCTGGGGGGTGTGCACGCGATAGCATTGCGAGACGCTGGA";
  std::string querySeq = "TCTATCACCCATCGTACCACTCACGGGAGCTCTCC";
  // As size of target sequence is less than query sequence,
  // required number of locations is zero.
  BOOST_REQUIRE_EQUAL(getQuerySeqMatchCount("AGCT", querySeq, 0.2), 0);
  // According to query sequence and target sequence, if it aligns from
  // location 9 of target sequence, then mismatch rate is ~0.14(<=0.2).
  // Mismatch rate in rest of location is more than 0.2.
  BOOST_REQUIRE_EQUAL(getQuerySeqMatchCount(targetSeq, querySeq, 0.2), 1);
  // According to query sequence and target sequence, if it aligns from
  // location 9 and 21 of target sequence, then mismatch rate is ~0.14(<=0.2)
  // and ~0.54(<=0.6) respectively. Mismatch rate in rest of location is more than 0.2.
  BOOST_REQUIRE_EQUAL(getQuerySeqMatchCount(targetSeq, querySeq, 0.6), 2);
}

// Test whether a candidate small SV alignment should be filtered due to low quality.
// Test the following cases:
// 1. minAlignRefSpan = 30 (simple SV) & 35 (complex SV)
// 2. minAlignReadLength = 30 (simple SV) & 35 (complex SV)
// 3. Fraction of alignment score should be greater than 0.75 where
//    fraction of alignment score is calculated as
//    alignment_score / (clipped_read_length * match score)
// 4. Soft clip should not be considered.
BOOST_AUTO_TEST_CASE(test_isLowQualitySmallSVAlignment)
{
  // Match score = 1, mismatch penalty = -4,
  // gap opening penalty = -6, gap extension penalty = -2,
  // clipping penalty = -5.
  AlignmentScores<int> scores(1, -4, -6, -2, -5);
  std::string          testCigar1("35=5I30=6D10=3I");
  ALIGNPATH::path_t    path1;
  cigar_to_apath(testCigar1.c_str(), path1);
  // According to penalty, alignment score = 29 and optimal score = 83.
  // So score fraction is 29/83 = ~0.34 (<0.75)
  BOOST_REQUIRE(isLowQualitySmallSVAlignment(100, scores, false, false, path1));

  std::string       testCigar2("35=5X30=6D10=3I");
  ALIGNPATH::path_t path2;
  cigar_to_apath(testCigar2.c_str(), path2);
  // According to penalty, alignment score = 25 and optimal score = 83.
  // So score fraction is 25/83 = ~0.30 (<0.75)
  BOOST_REQUIRE(isLowQualitySmallSVAlignment(100, scores, false, false, path2));

  std::string       testCigar3("35=");
  ALIGNPATH::path_t path3;
  cigar_to_apath(testCigar3.c_str(), path3);
  // According to penalty, alignment score = 35 and optimal score = 35.
  // So score fraction is 35/35 = 1.0 (>=0.75)
  BOOST_REQUIRE(!isLowQualitySmallSVAlignment(100, scores, false, false, path3));

  std::string       testCigar4("25=");
  ALIGNPATH::path_t path4;
  cigar_to_apath(testCigar4.c_str(), path4);
  // Here reference projection length is 25 which is less than minAlignRefSpan 30.
  BOOST_REQUIRE(isLowQualitySmallSVAlignment(100, scores, false, false, path4));

  std::string       testCigar5("25=23S");
  ALIGNPATH::path_t path5;
  cigar_to_apath(testCigar5.c_str(), path5);
  // Here clipped path size 25 which is less than minAlignReadLength 30.
  BOOST_REQUIRE(isLowQualitySmallSVAlignment(100, scores, false, false, path5));

  std::string       testCigar6("30=23S");
  ALIGNPATH::path_t path6;
  cigar_to_apath(testCigar6.c_str(), path6);
  // Here reference projection length and clipped path size satisfy the thresholds.
  // Soft-clip will not be considered in score calculation.
  // According to penalty, alignment score = 30 and optimal score = 30.
  // So score fraction is 30/30 = 1.0 (>=0.75)
  BOOST_REQUIRE(!isLowQualitySmallSVAlignment(100, scores, false, false, path6));

  // Following two cases are for complex candidate
  std::string       testCigar7("30=");
  ALIGNPATH::path_t path7;
  cigar_to_apath(testCigar7.c_str(), path4);
  // Here reference projection length is 30 which is less than minAlignRefSpan 35.
  BOOST_REQUIRE(isLowQualitySmallSVAlignment(100, scores, false, true, path7));

  std::string       testCigar8("30=23S");
  ALIGNPATH::path_t path8;
  cigar_to_apath(testCigar8.c_str(), path8);
  // Here clipped path size 30 which is less than minAlignReadLength 35.
  BOOST_REQUIRE(isLowQualitySmallSVAlignment(100, scores, false, true, path8));
}

// Test the following cases:
// 1. Whether a breakend insertion/deletion is inserted to SV candidate or not.
// 2. Insertion size = size of the candidate insert sequence
// 3. Deletion size = Difference between the two breakpoint start position - 1
// 4. All the above points are valid for indel type sv, otherwise it will not do anything.
// Following two schematic diagrams are vaild for indel type SV
//                   BP1               BP2
// Category-1   >>>>>>> -------------- <<<<<<< (BP1 Right-Open and BP2 Left-Open )
//
//                   BP2               BP1
// Category-2   >>>>>>> -------------- <<<<<<< (BP1 Left-Open and BP2 Right-Open )
BOOST_AUTO_TEST_CASE(test_addCigarToSpanningAlignment)
{
  // Category-1 (See docs)
  SVCandidate candidate1;
  candidate1.bp1.state    = SVBreakendState::RIGHT_OPEN;
  candidate1.bp2.state    = SVBreakendState::LEFT_OPEN;
  candidate1.bp1.interval = GenomeInterval(0, 1000, 1200);
  candidate1.bp2.interval = GenomeInterval(0, 2000, 2100);
  candidate1.insertSeq    = "AGAGCTAGT";
  // Insertion szie = 9
  // Deletion size = 2000 - 1000 - 1 = 999
  addCigarToSpanningAlignment(candidate1);
  BOOST_REQUIRE_EQUAL(candidate1.insertAlignment.size(), 2);
  BOOST_REQUIRE_EQUAL(candidate1.insertAlignment[0].type, ALIGNPATH::INSERT);
  BOOST_REQUIRE_EQUAL(candidate1.insertAlignment[0].length, 9);
  BOOST_REQUIRE_EQUAL(candidate1.insertAlignment[1].type, ALIGNPATH::DELETE);
  BOOST_REQUIRE_EQUAL(candidate1.insertAlignment[1].length, 999);

  // Category-2 (See docs)
  SVCandidate candidate2;
  candidate2.bp1.interval = GenomeInterval(0, 2000, 2100);
  candidate2.bp2.interval = GenomeInterval(0, 1000, 1200);
  candidate2.bp1.state    = SVBreakendState::LEFT_OPEN;
  candidate2.bp2.state    = SVBreakendState::RIGHT_OPEN;
  candidate2.insertSeq    = "AGAGCTAGT";
  // Insertion szie = 9
  // Deletion size = 2000 - 1000 - 1 = 999
  addCigarToSpanningAlignment(candidate2);
  BOOST_REQUIRE_EQUAL(candidate2.insertAlignment.size(), 2);
  BOOST_REQUIRE_EQUAL(candidate2.insertAlignment[0].type, ALIGNPATH::INSERT);
  BOOST_REQUIRE_EQUAL(candidate2.insertAlignment[0].length, 9);
  BOOST_REQUIRE_EQUAL(candidate2.insertAlignment[1].type, ALIGNPATH::DELETE);
  BOOST_REQUIRE_EQUAL(candidate2.insertAlignment[1].length, 999);

  // SV type is not indel. So it will not add any insertion
  // and deletion
  SVCandidate candidate3;
  candidate3.bp1.state    = SVBreakendState::RIGHT_OPEN;
  candidate3.bp2.state    = SVBreakendState::RIGHT_OPEN;
  candidate3.bp1.interval = GenomeInterval(0, 1000, 1200);
  candidate3.bp2.interval = GenomeInterval(0, 2000, 2100);
  candidate3.insertSeq    = "AGAGCTAGT";
  addCigarToSpanningAlignment(candidate3);
  BOOST_REQUIRE_EQUAL(candidate3.insertAlignment.size(), 0);
}

// Test whether max score contiguous path(starting from the beginning) of alignment
// satisfies the following criteria:
// 1. Contig length should be at least 40 in the max score contiguous path.
// 2. Reference projection length should be at least 40 in the max score contiguous path.
// 3. Fraction of alignment score relative to the optimal score in the max score path
//    should be greater than 0.75 where this alignment score fraction is calculated as
//    alignment_score / (clipped_contig_length * match score).
// 4. Length of the unaligned portion of contig means (read length -contig length) should be greater than 40.
//    For example:  If a cigar is 35=5I30=6D10=3I, then max score contiguous path
//    is 35=5I30=.
BOOST_AUTO_TEST_CASE(test_isLargeInsertSegment)
{
  // Match score = 1, mismatch penalty = -4,
  // gap opening penalty = -6, gap extension penalty = -2,
  // clipping penalty = -5.
  AlignmentScores<int> scores(1, -4, -6, -2, -5);
  AlignerBase<int>     alignerBase(scores);

  // contig length in the max score contiguous path
  unsigned contigOffset1;
  // Reference projection length in the max score contiguous path
  unsigned refoffset1;
  // alignment score in the max score contiguous path
  int               alignmentScore1;
  std::string       testCigar1("35=5I30=6D10=3I");
  ALIGNPATH::path_t path1;
  cigar_to_apath(testCigar1.c_str(), path1);
  // According to 35=5I30=6D10=3I, max score contiguous path is 35=5I30=
  // But here length of the unaligned portion of contig = 83 - 70 = 13
  // which is less than 40. So the api will return false.
  BOOST_REQUIRE(!isLargeInsertSegment(alignerBase, path1, contigOffset1, refoffset1, alignmentScore1));
  // Length of the offsets according to 35=5I30=
  BOOST_REQUIRE_EQUAL(contigOffset1, 70);
  BOOST_REQUIRE_EQUAL(refoffset1, 65);
  BOOST_REQUIRE_EQUAL(alignmentScore1, 49);

  unsigned          contigOffset2;
  unsigned          refoffset2;
  int               alignmentScore2;
  std::string       testCigar2("100=");
  ALIGNPATH::path_t path2;
  cigar_to_apath(testCigar2.c_str(), path2);
  // According to 100=, max score contiguous path is 100=
  // But here length of the unaligned portion of contig = 100 - 100 = 0
  // which is less than 40. So the api will return false.
  BOOST_REQUIRE(!isLargeInsertSegment(alignerBase, path2, contigOffset2, refoffset2, alignmentScore2));
  // Length of the offsets according to 100=
  BOOST_REQUIRE_EQUAL(contigOffset2, 100);
  BOOST_REQUIRE_EQUAL(refoffset2, 100);
  BOOST_REQUIRE_EQUAL(alignmentScore2, 100);

  unsigned          contigOffset3;
  unsigned          refoffset3;
  int               alignmentScore3;
  std::string       testCigar3("200=40X");  // maxpath is 200=
  ALIGNPATH::path_t path3;
  cigar_to_apath(testCigar3.c_str(), path3);
  // According to 200=40X, max score contiguous path is 200=
  // Here length of the unaligned portion of contig = 240 - 200 = 40 which is fine.
  // And the fraction of alignment score = 200 / 200 = 1 (>0.75).
  BOOST_REQUIRE(isLargeInsertSegment(alignerBase, path3, contigOffset3, refoffset3, alignmentScore3));
  // Length of the offsets according to 200=40X
  BOOST_REQUIRE_EQUAL(contigOffset3, 200);
  BOOST_REQUIRE_EQUAL(refoffset3, 200);
  BOOST_REQUIRE_EQUAL(alignmentScore3, 200);

  unsigned          contigOffset4;
  unsigned          refoffset4;
  int               alignmentScore4;
  std::string       testCigar4("49X200=40X");
  ALIGNPATH::path_t path4;
  cigar_to_apath(testCigar4.c_str(), path4);
  // According to 49X200=40X, max score contiguous path is 49X200=
  // Here length of the unaligned portion of contig = 289- 249 = 40 which is fine.
  // And the fraction of alignment score = 4 / 289 = ~0.01 (<0.75).
  BOOST_REQUIRE(!isLargeInsertSegment(alignerBase, path4, contigOffset4, refoffset4, alignmentScore4));
  // Length of the offsets according to 200=40X
  BOOST_REQUIRE_EQUAL(contigOffset4, 249);
  BOOST_REQUIRE_EQUAL(refoffset4, 249);
  BOOST_REQUIRE_EQUAL(alignmentScore4, 4);

  std::string       testCigar5("25=1I14=");
  ALIGNPATH::path_t path5;
  cigar_to_apath(testCigar5.c_str(), path5);
  // According to 25=1I14=, max score contiguous path is 25=1I14=
  // But here reference projection length is 39 which is less than 40.
  // So the api will return false.
  BOOST_REQUIRE(!isLargeInsertSegment(alignerBase, path5, contigOffset3, refoffset3, alignmentScore3));

  std::string       testCigar6("25=1D14=");
  ALIGNPATH::path_t path6;
  cigar_to_apath(testCigar6.c_str(), path6);
  // According to 25=1D14=, max score contiguous path is 25=1D14=
  // But here contig length is 39 which is less than 40.
  // So the api will return false.
  BOOST_REQUIRE(!isLargeInsertSegment(alignerBase, path6, contigOffset3, refoffset3, alignmentScore3));
}

// Test whether an alignment satisfies the test_isLargeInsertSegment
// either from the beginning or from the end or both.
BOOST_AUTO_TEST_CASE(test_isLargeInsertAlignment)
{
  // Match score = 1, mismatch penalty = -4,
  // gap opening penalty = -6, gap extension penalty = -2,
  // clipping penalty = -5.
  AlignmentScores<int> scores(1, -4, -6, -2, -5);
  GlobalAligner<int>   globalAligner(scores);
  LargeInsertionInfo   largeInsertionInfo1;
  std::string          testCigar1("200=40X");
  ALIGNPATH::path_t    path1;
  cigar_to_apath(testCigar1.c_str(), path1);
  // From start it is satisfied test_isLargeInsertSegment but from the end
  // it is not satisfied test_isLargeInsertSegment.
  BOOST_REQUIRE(isLargeInsertAlignment(globalAligner, path1, largeInsertionInfo1));
  BOOST_REQUIRE(largeInsertionInfo1.isLeftCandidate);
  BOOST_REQUIRE(!largeInsertionInfo1.isRightCandidate);

  LargeInsertionInfo largeInsertionInfo2;
  std::string        testCigar2("40X200=");
  ALIGNPATH::path_t  path2;
  cigar_to_apath(testCigar2.c_str(), path2);
  // From start it is not satisfied test_isLargeInsertSegment but from the end
  // it is satisfied test_isLargeInsertSegment.
  BOOST_REQUIRE(isLargeInsertAlignment(globalAligner, path2, largeInsertionInfo2));
  BOOST_REQUIRE(!largeInsertionInfo2.isLeftCandidate);
  BOOST_REQUIRE(largeInsertionInfo2.isRightCandidate);

  LargeInsertionInfo largeInsertionInfo3;
  std::string        testCigar3("150=");  // maxpath is 200=
  ALIGNPATH::path_t  path3;
  cigar_to_apath(testCigar3.c_str(), path3);
  // From start and end it is not satisfied test_isLargeInsertSegment.
  BOOST_REQUIRE(!isLargeInsertAlignment(globalAligner, path3, largeInsertionInfo3));
  BOOST_REQUIRE(!largeInsertionInfo3.isLeftCandidate);
  BOOST_REQUIRE(!largeInsertionInfo3.isRightCandidate);
}

// Test whether an alignment segment satisfies test_isLargeInsertSegment.
BOOST_AUTO_TEST_CASE(test_isFinishedLargeInsertAlignment)
{
  AlignmentScores<int> scores(1, -4, -6, -1, -5);
  GlobalAligner<int>   globalAligner(scores);
  std::string          testCigar1("40I350=40I");
  ALIGNPATH::path_t    path1;
  cigar_to_apath(testCigar1.c_str(), path1);
  BOOST_REQUIRE(isFinishedLargeInsertAlignment(globalAligner, path1, std::pair<unsigned, unsigned>(0, 2), 0));
  BOOST_REQUIRE(
      !isFinishedLargeInsertAlignment(globalAligner, path1, std::pair<unsigned, unsigned>(0, 2), 1));
}

// Test the original location of masked position where masked position
// is constructed using exclusion_block data structure which conatins:
// 1. start of the excluded region
// 2. Number of Bp excluded(l).
// 3. Number of 'N' added in place of excluded sequence(N).
// Consider the following schematic diagram of reduced coordinate:
// | l = 10, N=3 |    | l = 10, N=3 |    | l = 10, N=3 |
// --------------------------------------------------------------
// | <- excl-1 ->|  P | <- excl-2 ->|    | <- excl-3 ->|
// We need to find the actual coordinate of P. So for every block need
// to add offset (l-N) to P and check whether it is greater than the start
// of the next block and so on.
BOOST_AUTO_TEST_CASE(test_translateMaskedPos)
{
  // start of excluded region = 1
  // Length of the excluded region = 10
  // Number of 'N' added in place of excluded sequence = 3
  exclusion_block block1(1, 10, 3);
  // start of excluded region = 12
  // Length of the excluded region = 10
  // Number of 'N' added in place of excluded sequence = 3
  exclusion_block block2(12, 10, 3);
  // start of excluded region = 22
  // Length of the excluded region = 10
  // Number of 'N' added in place of excluded sequence = 3
  exclusion_block              block3(22, 10, 3);
  std::vector<exclusion_block> blocks = {block1, block2, block3};
  // Original coordinate = 5 + (10-3) + (10-3) = 19 (Based on
  // the above diagram)
  BOOST_REQUIRE_EQUAL(translateMaskedPos(blocks, 5), 19);
}

// Test the breakend interval range to be reported in the output vcf.
// Test the following cases:
// 1. In case of sv breakend state is RIGHT_OPEN
// 2. In case of sv breaked state is LEFT_OPEN
BOOST_AUTO_TEST_CASE(test_adjustAssembledBreakend)
{
  reference_contig_segment referenceContigSegment;
  referenceContigSegment.seq() =
      "GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGT"
      "ATTTTCGTCTGGGGGGTGTGCACGCGATAGCATTGCGAGACGCTGGA";
  Alignment alignment;
  alignment.beginPos = 9;
  std::string testCigar("35M");
  cigar_to_apath(testCigar.c_str(), alignment.apath);
  SVBreakend svBreakend;
  svBreakend.state = SVBreakendState::RIGHT_OPEN;
  // homologous jump range across the breakend = 2
  // As state is RIGHT_OPEN, range start is alignment end -1 = 44 - 1 = 43
  // and alignment end = range start + jump range + 1 = 46
  adjustAssembledBreakend(alignment, true, 2, referenceContigSegment, false, svBreakend);
  BOOST_REQUIRE_EQUAL(svBreakend.interval.range.begin_pos(), 43);
  BOOST_REQUIRE_EQUAL(svBreakend.interval.range.end_pos(), 46);

  svBreakend.state = SVBreakendState::LEFT_OPEN;
  // homologous jump range across the breakend = 2
  // As state is LEFT_OPEN, range start is alignment start -jump range = 9 - 2 = 7
  // and alignment end = alignment start + 1 = 10
  adjustAssembledBreakend(alignment, true, 2, referenceContigSegment, false, svBreakend);
  BOOST_REQUIRE_EQUAL(svBreakend.interval.range.begin_pos(), 7);
  BOOST_REQUIRE_EQUAL(svBreakend.interval.range.end_pos(), 10);
}

// Test whether a candidate spanning SV alignment should be filtered due to low quality.
// Test the following cases:
// 1. minAlignReadLength = 30 (DNA) & 20 (RNA)
// 2. Fraction of alignment score relative to the optimal score should be greater than 0.75 where
//    the fraction of alignment score is calculated as
//    alignment_score / (clipped_read_length * match score)
// 3. Soft clip should not be considered.
BOOST_AUTO_TEST_CASE(test_isLowQualitySpanningSVAlignment)
{
  // Match score = 1, mismatch penalty = -4,
  // gap opening penalty = -6, gap extension penalty = -2,
  // clipping penalty = -5.
  AlignmentScores<int> scores(1, -4, -6, -2, -5);
  std::string          testCigar1("35=5I30=6D10=3I");  // match, insertion and deletion
  ALIGNPATH::path_t    path1;
  cigar_to_apath(testCigar1.c_str(), path1);
  // According to penalty, alignment score = 29 and optimal score = 83.
  // So score fraction is 29/83 = ~0.34 (<0.75)
  BOOST_REQUIRE(isLowQualitySpanningSVAlignment(100, scores, false, false, path1));

  std::string       testCigar2("35=5X30=6D10=3I");  // match, mismatch, insertion and deletion
  ALIGNPATH::path_t path2;
  cigar_to_apath(testCigar2.c_str(), path2);
  // According to penalty, alignment score = 29 and optimal score = 83.
  // So score fraction is 25/83 = ~0.30 (<0.75)
  BOOST_REQUIRE(isLowQualitySpanningSVAlignment(100, scores, false, false, path2));

  std::string       testCigar3("35=");
  ALIGNPATH::path_t path3;
  cigar_to_apath(testCigar3.c_str(), path3);
  // According to penalty, alignment score = 35 and optimal score = 35.
  // So score fraction is 35/35 = 1.0 (>=0.75)
  BOOST_REQUIRE(!isLowQualitySpanningSVAlignment(100, scores, false, false, path3));

  std::string       testCigar4("25=");
  ALIGNPATH::path_t path4;
  cigar_to_apath(testCigar4.c_str(), path4);
  // Here read size 25 which is less than minAlignReadLength 30.
  BOOST_REQUIRE(isLowQualitySpanningSVAlignment(100, scores, false, false, path4));

  std::string       testCigar5("25=23S");
  ALIGNPATH::path_t path5;
  cigar_to_apath(testCigar5.c_str(), path5);
  // Here clipped path size 25 which is less than minAlignReadLength 30.
  BOOST_REQUIRE(isLowQualitySpanningSVAlignment(100, scores, false, false, path5));

  std::string       testCigar6("30=23S");
  ALIGNPATH::path_t path6;
  cigar_to_apath(testCigar6.c_str(), path6);
  // According to penalty, non-clipped alignment score = 30 and optimal score = 30.
  // So score fraction is 30/30 = 1.0 (>=0.75)
  BOOST_REQUIRE(!isLowQualitySpanningSVAlignment(100, scores, false, false, path6));

  // Following two cases are for RNA
  std::string       testCigar7("15=");
  ALIGNPATH::path_t path7;
  cigar_to_apath(testCigar7.c_str(), path4);
  // Here read size 15 which is less than minAlignReadLength 20.
  BOOST_REQUIRE(isLowQualitySpanningSVAlignment(100, scores, false, true, path7));

  std::string       testCigar8("15=23S");
  ALIGNPATH::path_t path8;
  cigar_to_apath(testCigar8.c_str(), path8);
  // Here clipped path size 15 which is less than minAlignReadLength 20.
  BOOST_REQUIRE(isLowQualitySpanningSVAlignment(100, scores, false, true, path8));

  std::string       testCigar9("25=23S");
  ALIGNPATH::path_t path9;
  cigar_to_apath(testCigar9.c_str(), path9);
  // According to penalty, non-clipped alignment score = 25 and optimal score = 25.
  // So score fraction is 25/25 = 1.0 (>=0.75)
  BOOST_REQUIRE(!isLowQualitySpanningSVAlignment(100, scores, false, true, path9));
}

// After getting the jump alignment (using assembly) for a spanning candidate SV,
// adjust the breakend region of SV candidate based on the original
// genome coordinate. Test the breakend interval range to be reported in the output vcf.
// Test the following cases:
// 1. In case of sv breakend state is RIGHT_OPEN
// 2. In case of sv breaked state is LEFT_OPEN
BOOST_AUTO_TEST_CASE(test_generateRefinedSVCandidateFromJumpAlignment)
{
  // Temporary assembly data
  SVCandidateAssemblyData svCandidateAssemblyData;
  // index of the alignment
  svCandidateAssemblyData.bestAlignmentIndex = 0;
  // Created two alignments with jump size = 2
  Alignment alignment1;
  alignment1.beginPos = 40;
  std::string testCigar1("35=");
  cigar_to_apath(testCigar1.c_str(), alignment1.apath);
  Alignment alignment2;
  alignment2.beginPos = 70;
  std::string testCigar2("10=75I1D2=");
  cigar_to_apath(testCigar2.c_str(), alignment2.apath);
  SVCandidateAssemblyData::JumpAlignmentResultType jumpAlignmentResultType;
  jumpAlignmentResultType.align1    = alignment1;
  jumpAlignmentResultType.align2    = alignment2;
  jumpAlignmentResultType.jumpRange = 2;
  svCandidateAssemblyData.spanningAlignments.push_back(jumpAlignmentResultType);
  svCandidateAssemblyData.bp1ref.seq() =
      "TCTATCACCCATCGTACCACTCACGGGAGCTCTCCTCTATCACCCATCGTACCACTCACGGGAGCTCTCC"
      "TCTATCACCCATCGTACCACTCACGGGAGCTCTCCTCTATCACCCATCGTACCACTCACGGGAGCTCTCC"
      "TCTATCACCCATCGTACCACTCACGGGAGCTCTCC";
  svCandidateAssemblyData.bp2ref.seq() =
      "GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGT"
      "ATTTTCGTCTGGGGGGTGTGCACGCGATAGCATTGCGAGACGCTGGA"
      "GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGT"
      "ATTTTCGTCTGGGGGGTGTGCACGCGATAGCATTGCGAGACGCTGGA"
      "GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGT"
      "ATTTTCGTCTGGGGGGTGTGCACGCGATAGCATTGCGAGACGCTGGA";

  svCandidateAssemblyData.bp1ref.set_offset(100);
  svCandidateAssemblyData.bp2ref.set_offset(150);
  // BP1 and BP2 both are RIGHT_OPEN
  // As state is RIGHT_OPEN,
  //     BP1:
  //         Range start = bp1_reffset+ alignment end - 1 = 175 - 1 = 174
  //         Rand end = range start + jump range + 1 = 177
  //     BP2 (Reverse strand):
  //         Range start = bp2_refoffset + (bp2_ref_size - alignment start - 1) - jump size
  //                     = 150 + (306 - 70 - 1) - 2 = 383
  //         Rand end = bp2_refoffset + (bp2_ref_size - alignment start - 1) + 1
  //                  = 150 + (306 - 70 - 1) + 1 = 386
  SVCandidate svCandidate1;
  svCandidate1.bp1.state = SVBreakendState::RIGHT_OPEN;
  svCandidate1.bp2.state = SVBreakendState::RIGHT_OPEN;
  // As here BP2 is right opened that means BP2 is reversed.
  svCandidateAssemblyData.bporient.isBp2Reversed = true;
  generateRefinedSVCandidateFromJumpAlignment(svCandidateAssemblyData, svCandidate1);
  // results should be according to test_adjustAssembledBreakend
  BOOST_REQUIRE_EQUAL(svCandidate1.bp1.interval.range.begin_pos(), 174);
  BOOST_REQUIRE_EQUAL(svCandidate1.bp1.interval.range.end_pos(), 177);
  BOOST_REQUIRE_EQUAL(svCandidate1.bp2.interval.range.begin_pos(), 383);
  BOOST_REQUIRE_EQUAL(svCandidate1.bp2.interval.range.end_pos(), 386);

  // Created two alignments with jump size = 0
  alignment1.beginPos = 397;
  std::string testCigar3("92=");
  cigar_to_apath(testCigar3.c_str(), alignment1.apath);
  alignment2.beginPos = 510;
  std::string testCigar4("110=75I1D2=");
  cigar_to_apath(testCigar4.c_str(), alignment2.apath);
  jumpAlignmentResultType.align1 = alignment1;
  jumpAlignmentResultType.align2 = alignment2;
  // jump range = 2
  jumpAlignmentResultType.jumpRange = 2;
  svCandidateAssemblyData.spanningAlignments.clear();
  svCandidateAssemblyData.spanningAlignments.push_back(jumpAlignmentResultType);
  // As BP2 is LEFT_OPEN, BP2 is not in reverse strand.
  svCandidateAssemblyData.bporient.isBp2Reversed = false;
  // BP1 is RIGHT_OPEN and BP2 is LEFT_OPEN
  //     BP1:
  //         Range start = bp1_refoffset + alignment end - 1 = 589 - 1 = 588
  //         Rand end = range start + jump range + 1 = 591
  //     BP2:
  //         Range start = bp2_refoffset + alignment start = 150 + 510 = 660
  //         Rand end = range start  + jump range + 1 = 663
  SVCandidate svCandidate2;
  svCandidate2.bp1.state = SVBreakendState::RIGHT_OPEN;
  svCandidate2.bp2.state = SVBreakendState::LEFT_OPEN;
  generateRefinedSVCandidateFromJumpAlignment(svCandidateAssemblyData, svCandidate2);
  // results should be according to test_adjustAssembledBreakend
  BOOST_REQUIRE_EQUAL(svCandidate2.bp1.interval.range.begin_pos(), 588);
  BOOST_REQUIRE_EQUAL(svCandidate2.bp1.interval.range.end_pos(), 591);
  BOOST_REQUIRE_EQUAL(svCandidate2.bp2.interval.range.begin_pos(), 660);
  BOOST_REQUIRE_EQUAL(svCandidate2.bp2.interval.range.end_pos(), 663);
}

// Test whether jump alignment QC fails or not.
// Test the following cases:
// 1. Cigar should be non-empty
// 2. Minimum read length should be 20
// 3. QC should fail if one of the alignment fails one of the above two condidtions.
BOOST_AUTO_TEST_CASE(test_isJumpAlignmentQCFail)
{
  Alignment alignment1;
  alignment1.beginPos = 406;
  std::string testCigar1("94=");
  cigar_to_apath(testCigar1.c_str(), alignment1.apath);
  Alignment alignment2;
  alignment2.beginPos = 510;
  std::string testCigar2("110=75I1D2=");
  cigar_to_apath(testCigar2.c_str(), alignment2.apath);
  SVCandidateAssemblyData::JumpAlignmentResultType jumpAlignmentResultType;
  jumpAlignmentResultType.align1 = alignment1;
  jumpAlignmentResultType.align2 = alignment2;
  // Both alignment1 and alignment2 satisfy point-1 and point-2
  BOOST_REQUIRE(!isJumpAlignmentQCFail(jumpAlignmentResultType));

  std::string testCigar3("");
  jumpAlignmentResultType.align1.beginPos = -1;
  cigar_to_apath(testCigar3.c_str(), jumpAlignmentResultType.align1.apath);
  // Alignment-2 satisfies point-1 and point-2 but alignment-1 doest not satisfy point-1
  BOOST_REQUIRE(isJumpAlignmentQCFail(jumpAlignmentResultType));

  std::string testCigar4("19=");
  jumpAlignmentResultType.align1.beginPos = 406;
  cigar_to_apath(testCigar4.c_str(), jumpAlignmentResultType.align1.apath);
  // Alignment-2 satisfies point-1 and point-2 but alignment-1 doest not satisfy point-2
  BOOST_REQUIRE(isJumpAlignmentQCFail(jumpAlignmentResultType));

  std::string testCigar5("92=");
  jumpAlignmentResultType.align1.beginPos = 406;
  cigar_to_apath(testCigar5.c_str(), jumpAlignmentResultType.align1.apath);
  std::string testCigar6("");
  jumpAlignmentResultType.align2.beginPos = -1;
  cigar_to_apath(testCigar6.c_str(), jumpAlignmentResultType.align2.apath);
  // Alignment-1 satisfies point-1 and point-2 but alignment-2 doest not satisfy point-1
  BOOST_REQUIRE(isJumpAlignmentQCFail(jumpAlignmentResultType));

  std::string testCigar7("19=");
  jumpAlignmentResultType.align2.beginPos = 406;
  cigar_to_apath(testCigar7.c_str(), jumpAlignmentResultType.align2.apath);
  // Alignment-1 satisfies point-1 and point-2 but alignment-2 doest not satisfy point-2
  BOOST_REQUIRE(isJumpAlignmentQCFail(jumpAlignmentResultType));

  std::string testCigar8("");
  jumpAlignmentResultType.align1.beginPos = -1;
  cigar_to_apath(testCigar8.c_str(), jumpAlignmentResultType.align1.apath);
  std::string testCigar9("");
  jumpAlignmentResultType.align2.beginPos = -1;
  cigar_to_apath(testCigar9.c_str(), jumpAlignmentResultType.align2.apath);
  // Alignment-1 and alignment-2 do not satisfy point-1 and point-2
  BOOST_REQUIRE(isJumpAlignmentQCFail(jumpAlignmentResultType));
}

// The motivation behind this api is:
// Let's say indel segments in a cigar are (a1,b1), (a2,b2) and (a3,b3) and all these indel stretches are
// satisfied min indel threshold criteria. API tries to check whether the segment before any of these indel
// stretch is high quality or not. Similarly it tries to check whether the segment after
// any of these indel stretch is high quality or not. If the segment before or after any of these indel is low
// quality drop that indel from the indel segments. After dropping those invalid indel segments, if the
// segment still contains some indel segments which satisfy min indel threshold criteria it returns true. A
// segment is said to be high quality if the
// following conditions are satisfied:
// 1. minAlignRefSpan = 30 (simple SV) & 35 (complex SV)
// 2. minAlignReadLength = 30 (simple SV) & 35 (complex SV)
// 3. Fraction of alignment score should be greater than 0.75 where
//    fraction of alignment score is calculated as
//    alignment_score / (clipped_read_length * match score)

// Test the following cases for candidate generation from complex SV contig alignment:
// 1. If minimum candidate indel is not present in the read alignment, it will return false.
// 2. If number of larger indel candidate segments is 1 and it is Low Quality SmallSV Alignment (as explained
// in
//    test_ISLowQualitySmallSVAlignment), it will return false.
// 3. After discarding invalid indel segments(as described above) if it still contains some indel segments
// which are
//    more than min indel threshold, it will return true.
BOOST_AUTO_TEST_CASE(test_findCandidateVariantsFromComplexSVContigAlignment)
{
  std::string contigSeq =
      "TCTATCACCCATCGTACCACTCACGGGAGCTCTCCTCTATCACCCATCGTACCACTCACGGGAGCTCTCC"
      "TCTATCACCCATCGTACCACTCACGGGAGCTCTCCTCTATCACCCATCGTACCACTCACGGGAGCTCTCC"
      "TCTATCACCCATCGTACCACTCACGGGAGCTCTCC";
  std::string refSeq =
      "GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGT"
      "ATTTTCGTCTGGGGGGTGTGCACGCGATAGCATTGCGAGACGCTGGA"
      "GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGT"
      "ATTTTCGTCTGGGGGGTGTGCACGCGATAGCATTGCGAGACGCTGGA"
      "GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGT"
      "ATTTTCGTCTGGGGGGTGTGCACGCGATAGCATTGCGAGACGCTGGA";

  std::vector<std::pair<unsigned, unsigned>> candidateSegments;

  // Match score = 1, mismatch penalty = -4,
  // gap opening penalty = -6, gap extension penalty = -2,
  // clipping penalty = -5.
  AlignmentScores<int> scores(1, -4, -6, -2, -5);
  std::string          testCigar1("35=5I30=6D10=3I");  // match, insertion and deletion
  Alignment            alignment1;
  cigar_to_apath(testCigar1.c_str(), alignment1.apath);
  // Case-1 is designed here. Minimum indel size = 10. But cigar(35=5I30=6D10=3I) does not have any indel
  // with size more than 10.
  BOOST_REQUIRE(!findCandidateVariantsFromComplexSVContigAlignment(
      100, scores, alignment1, contigSeq, refSeq, 10, candidateSegments));

  // Case-2 is designed here. Min indel size = 6. So only one segment (6D) more than this.
  // So before this segment, cigar is 35=5I30= which is low quality as alignment score 0.75(49/70)
  // which is less than 0.75.
  BOOST_REQUIRE(!findCandidateVariantsFromComplexSVContigAlignment(
      100, scores, alignment1, contigSeq, refSeq, 6, candidateSegments));

  // Min indel threshold is 3. So Indel segments are [(0,0), (2,2), (4,4), (6,6)].
  // Before 1st 3I(0,0), there is no segment, so the segment before (0,0) is low quality.
  // Similarly the segment after last 3I(6,6) is low quality. So drop these two segments
  // [(0,0) & (6,6)] from indel segments.
  // After dropping both the segments from indel segments, still it contains 5I and 6D with
  // sizes larger than min indel threshold = 3.
  std::string testCigar2("3I100=5I30=6D100=3I");  // match, insertion and deletion
  Alignment   alignment2;
  cigar_to_apath(testCigar2.c_str(), alignment2.apath);
  alignment2.beginPos = 50;
  BOOST_REQUIRE(findCandidateVariantsFromComplexSVContigAlignment(
      100, scores, alignment2, contigSeq, refSeq, 3, candidateSegments));
}

// The api constructs a usable candidate SV from a smallIndel alignment section.
// It first identifies from alignment cigar a segment of interest with potential SV candidate,
// then set the breakpoint range by expanding it around the breakpoint without changing alignment
// score (tested by test_GetVariantRange). Test the expanded breakend region around bp1 and bp2.
BOOST_AUTO_TEST_CASE(test_setSmallCandSV)
{
  reference_contig_segment reference;
  reference.seq() =
      "AATATTACCCTTAGTTCTCATTAGTAAATGGATGTGAAAGTCAGGACATTCCAAAAGAAAGAACTAGAACTCACTCGGCC"
      "AGAAAACCCCCATTTCAGTTTTTACACAGAAAAATTCTTACAGTCTATGTTTCACTAAGAATGTCTGCTGTGCAAAACCCT"
      "CAAACTTTTTAGAACGTTTTTTTTGTTCCAAGTTAGAGAACGGCAATCAGTAATCTATTACCCAAAGTGCTTCTCCTTTCCA"
      "GGTTTCATGTTAGAGTGATTCTAATATGTGTGTGCTATCAACTGCCTACACAGAAAACTGAGAGACAAAGGCTTTCTCCTTTT"
      "CCACACATTATCCTTCATTCAGACTTAATGCCTGCAGGTCCGGTTTAATGATTTCCCAGAGTTTATGAACGAAAAAGAAAAAC"
      "AAAGCAATAAAAACAAAAATCAAAGTTAAATTTCCTCAAAAGTTTTCAAGAAGGAAGTAGTCAGGACAAAAACAAAGGGAATG"
      "AGGGCACTTTGTCTCAGGATACAAATTAAAGATCACTGTGGTGGCCTCTGTGGGGTGGTTATAAAGGGGACCAGGTGTATACT"
      "AGGAAGTCATTTAGTTTTAGAAATGTAAATATGTGTAAATGTTTTAATTTTACTCAACTTGCCAGAGGTAGAATGTCCCTGGA"
      "CAACTAACTGATACATTTCTTTAGGGCCAATCGCTGGCTTTAGAAGAGCCTCAGCTAATCACAGTAGAGCTGGACTGTTGTGG"
      "TTTTCCATTCCTTTGCATCGTATTCCTCAGTCTCTGCGGAAGGCACTGCTCCTTCCTTTCCTTTCTAAATCTCTCTCTGTCTCTC"
      "TCTCTCTCTCTCTCTCTCTCGCTCTCTCCTCCCCTAGTTTATCCTGGACTCATGCTGAGCTCAGCAACCCTTGAACTCATTTTCT"
      "ATCTGACATGT";
  const std::string readSequence(
      "ATTCCTTTGCATCGTATTCCTCAGTCTCTGCGGAAGGCACTGCTCCTTCCTTTCCTTTCTAAATCTCTCTCTGTCTCT"
      "CTCTCTCTCTCTCTCTCTCTCTCTCTCGCTCTCTCCTCCCCTAGTTTATCCTGGACTCATGCTGAGCTCAGCAACCCTT"
      "GAACTCATTTTCTATCTGAC");
  // If an alignment cigar is 73=6I98=, our segment interest is 6I, api constructs a sv region
  // around that 6I. It will discard 73= and 98= segments. Let's say reference range for 6I is [x, y).
  // While doing this computation at the corresponding reference range, api checks that how far it can
  // expand  around the breakpoint without changing alignment score. Let's say this range is [a, b)
  // So the SV breakend locations are:
  // bp1 = [x, x + b + 1)
  // bp2 = [y, y + b + 1)
  Alignment alignment;
  alignment.beginPos = 747;
  std::string cigar("73=6I98=");
  cigar_to_apath(cigar.c_str(), alignment.apath);
  SVCandidate candidate;
  GSCOptions  options;
  options.isOutputContig = true;
  std::pair<unsigned, unsigned> segment(std::pair<unsigned, unsigned>(1, 1));
  setSmallCandSV(reference, readSequence, alignment, segment, candidate, options);
  // ref region is [747 +73-1, 747 + 73) = [819, 820)
  // insertion location is [73, 74)
  // from 819 in reference coordinate and from 73 in read coordinate 26 read bases are matching with
  // reference bases. So final regions are [819, 819+26+1) = [819, 846) and
  // [820, 820+26+1) = [820, 847)
  BOOST_REQUIRE_EQUAL(candidate.bp1.interval, GenomeInterval(0, 819, 846));
  BOOST_REQUIRE_EQUAL(candidate.bp2.interval, GenomeInterval(0, 820, 847));
  BOOST_REQUIRE_EQUAL(candidate.insertSeq.size(), 6);  // for 6I
}

// Test the original location of masked position where masked position
// is constructed using exclusion_block data structure which conatins:
// 1. start of the excluded region
// 2. Number of Bp excluded(l).
// 3. Number of 'N' added in place of excluded sequence(N).
// Consider the following schematic diagram of reduced coordinate:
// | l = 10, N=3 |    | l = 10, N=3 |    | l = 10, N=3 |
// --------------------------------------------------------------
// | <- excl-1 ->|  P | <- excl-2 ->|    | <- excl-3 ->|
// We need to find the actual coordinate of P. So for every block need
// to add offset (l-N) to P and check whether it is greater than the start
// of the next block and so on.
BOOST_AUTO_TEST_CASE(test_translateMaskedAlignment)
{
  // start of excluded region = 1
  // Length of the excluded region = 10
  // Number of 'N' added in place of excluded sequence = 3
  exclusion_block block1(1, 10, 3);
  // start of excluded region = 22
  // Length of the excluded region = 10
  // Number of 'N' added in place of excluded sequence = 3
  exclusion_block block2(22, 10, 3);
  // start of excluded region = 34
  // Length of the excluded region = 10
  // Number of 'N' added in place of excluded sequence = 3
  exclusion_block              block3(34, 10, 3);
  std::vector<exclusion_block> blocks = {block1, block2, block3};

  std::string testCigar1("20D4=4D12=4I");  // match, insertion and deletion
  Alignment   alignment1;
  alignment1.beginPos = 2;
  cigar_to_apath(testCigar1.c_str(), alignment1.apath);
  // Based on the above schematic diagram this is the expected original cigar.
  std::string       expectedCIGAR("34D4=4D12=4I");
  ALIGNPATH::path_t expectedPath;
  cigar_to_apath(expectedCIGAR.c_str(), expectedPath);
  BOOST_REQUIRE(translateMaskedAlignment(alignment1, blocks));
  BOOST_REQUIRE_EQUAL(alignment1.apath, expectedPath);
}

// Test whether jump alignment result is low quality
// Jump alignment represents alignment of a query sequence which can switch over
// from reference1 to reference2. If either of the alignment is LowQualitySpanningSVAlignment
// then api will return true.
// isLowQualitySpanningSVAlignment is described in test_isLowQualitySpanningSVAlignment.
// Test the following cases:
// 1. isLowQualitySpanningSVAlignment is true for Alignment-1, not true for Alignment-2
// 2. isLowQualitySpanningSVAlignment is true for Alignment-2, not true for Alignment-1
// 3. isLowQualitySpanningSVAlignment is true for both Alignment-1 and Alignment-2
// 4. isLowQualitySpanningSVAlignment is not true for both Alignment-1 and Alignment-2
BOOST_AUTO_TEST_CASE(test_isLowQualityJumpAlignment)
{
  JumpAlignmentResult<int> alignmentResult;
  // Match score = 1, mismatch penalty = -4,
  // gap opening penalty = -6, gap extension penalty = -2,
  // clipping penalty = -5.
  AlignmentScores<int> scores(1, -4, -6, -2, -5);
  std::string          testCigar1("35=5I30=6D10=3I");  // match, insertion and deletion
  cigar_to_apath(testCigar1.c_str(), alignmentResult.align1.apath);
  std::string testCigar2("35=");
  cigar_to_apath(testCigar2.c_str(), alignmentResult.align2.apath);
  // Case-1 is designed here Where fraction of alignment score relative
  // to optimal score of alignment-1 is 0.349398 (<0.75) but fraction of
  // alignment score of alignment-2 is 1(>0.75).
  // So here alignment-1 is low quality and alignment-2 is high quality.
  BOOST_REQUIRE(isLowQualityJumpAlignment(alignmentResult, scores, false));
  // Case-2 is designed here Where fraction of alignment score relative
  // to optimal score of alignment-2 is 0.349398 (<0.75) but fraction of
  // alignment score of alignment-1 is 1(>0.75).
  // So here alignment-2 is low quality and alignment-1 is high quality.
  cigar_to_apath(testCigar2.c_str(), alignmentResult.align1.apath);
  cigar_to_apath(testCigar1.c_str(), alignmentResult.align2.apath);
  BOOST_REQUIRE(isLowQualityJumpAlignment(alignmentResult, scores, false));

  // Case-3 is designed here Where fraction of alignment score relative
  // to optimal score of alignment-1 is 0.349398 (<0.75) but fraction of
  // alignment score of alignment-2 is 0.349398(<0.75).
  // So here alignment-1 and alignment-2 both are low quality.
  cigar_to_apath(testCigar1.c_str(), alignmentResult.align1.apath);
  cigar_to_apath(testCigar1.c_str(), alignmentResult.align2.apath);
  BOOST_REQUIRE(isLowQualityJumpAlignment(alignmentResult, scores, false));

  // Case-4 is designed here Where fraction of alignment score relative
  // to optimal score of alignment-1 is 1 (>0.75) but fraction of
  // alignment score of alignment-2 is 1(>0.75).
  // So here alignment-1 and alignment-2 both are high quality.
  cigar_to_apath(testCigar2.c_str(), alignmentResult.align1.apath);
  cigar_to_apath(testCigar2.c_str(), alignmentResult.align2.apath);
  BOOST_REQUIRE(!isLowQualityJumpAlignment(alignmentResult, scores, false));
}

// Test the following cases:
// 1. Filter breakpoint contigs for large SV candidates as explained in test_isLowQualityJumpAlignment
// 2. Select the 'best' one based on alignment score and check alignment on selected contig
BOOST_AUTO_TEST_CASE(test_selectJumpContigDNA)
{
  // Match score = 1, mismatch penalty = -4,
  // gap opening penalty = -6, gap extension penalty = -2,
  // clipping penalty = -5.
  AlignmentScores<int>     scores(1, -4, -6, -2, -5);
  JumpAlignmentResult<int> alignmentResult1;
  std::string              testCigar1("35=5I30=6D10=3I");             // match, insertion and deletion
  cigar_to_apath(testCigar1.c_str(), alignmentResult1.align1.apath);  // alignment score is 29
  std::string testCigar2("35=");
  cigar_to_apath(testCigar2.c_str(), alignmentResult1.align2.apath);  // alignment score is 35
  alignmentResult1.score = 64;

  JumpAlignmentResult<int> alignmentResult2;
  std::string              testCigar3("35=5X30=6D10=3I");
  cigar_to_apath(testCigar3.c_str(), alignmentResult2.align1.apath);  // alignment score is 25
  std::string testCigar4("35=");
  cigar_to_apath(testCigar4.c_str(), alignmentResult2.align2.apath);  // alignment score is 35
  alignmentResult1.score = 60;

  // Case-1 is designed here. Out of two jump alignments, alignmentResult1 has max alignment score
  // But it is LowQualityJumpAlignment as one of the alignments (35=5I30=6D10=3I) is low quality.
  // The reson behind low quality is its alignment score is 0.349398 (29/83) which is less than 0.75.
  SVCandidateAssemblyData candidateAssemblyData;
  candidateAssemblyData.contigs.resize(2);
  candidateAssemblyData.spanningAlignments.push_back(alignmentResult1);
  candidateAssemblyData.spanningAlignments.push_back(alignmentResult2);
  BOOST_REQUIRE(!selectJumpContigDNA(candidateAssemblyData, scores));

  // Case-2 is designed here. Out of three jump alignments, alignmentResult3 has max alignment score.
  // And also alignmentResult3 has high quality jump alignment bacause all the
  // jump alignments' score are 1 (>0.75).
  JumpAlignmentResult<int> alignmentResult3;
  std::string              testCigar5("35=");
  cigar_to_apath(testCigar5.c_str(), alignmentResult3.align1.apath);  // alignment score is 35
  std::string testCigar6("35=");
  cigar_to_apath(testCigar6.c_str(), alignmentResult3.align2.apath);  // alignment score is 35
  alignmentResult3.score = 70;
  candidateAssemblyData.contigs.resize(3);
  candidateAssemblyData.spanningAlignments.push_back(alignmentResult3);
  BOOST_REQUIRE(selectJumpContigDNA(candidateAssemblyData, scores));
  BOOST_REQUIRE_EQUAL(candidateAssemblyData.bestAlignmentIndex, 2);
}

// Test the following cases for selection of contig in RNA:
// 1. Filter breakpoint contigs for large SV candidates as explained in test_isLowQualityJumpAlignment.
// 2. Select the 'best' one based on alignment score and check alignment on selected contig.
// 3. If two or more jump alignments support high quality spanning SV, best alignment index will be that
//    jump alignment whose number of supporting reads for contig is maximum.
BOOST_AUTO_TEST_CASE(test_selectJumpContigRNA)
{
  // Match score = 1, mismatch penalty = -4,
  // gap opening penalty = -6, gap extension penalty = -2,
  // clipping penalty = -5.
  AlignmentScores<int>     scores(1, -4, -6, -2, -5);
  JumpAlignmentResult<int> alignmentResult1;
  std::string              testCigar1("35=5I30=6D10=3I");             // match, insertion and deletion
  cigar_to_apath(testCigar1.c_str(), alignmentResult1.align1.apath);  // alignment score is 29
  std::string testCigar2("35=");
  cigar_to_apath(testCigar2.c_str(), alignmentResult1.align2.apath);  // alignment score is 35
  alignmentResult1.score = 64;

  JumpAlignmentResult<int> alignmentResult2;
  std::string              testCigar3("35=5X30=6D10=3I");
  cigar_to_apath(testCigar3.c_str(), alignmentResult2.align1.apath);  // alignment score is 25
  std::string testCigar4("35=");
  cigar_to_apath(testCigar4.c_str(), alignmentResult2.align2.apath);  // alignment score is 35
  alignmentResult1.score = 60;

  // Case-1 is designed here. Out of two jump alignments, alignmentResult1 has max alignment score
  // But it is LowQualityJumpAlignment as one of the alignments (35=5I30=6D10=3I) is low quality.
  // The reson behind low quality is its alignment score is 0.349398 (29/83) which is less than 0.75.
  SVCandidateAssemblyData candidateAssemblyData;
  candidateAssemblyData.contigs.resize(2);
  candidateAssemblyData.spanningAlignments.push_back(alignmentResult1);
  candidateAssemblyData.spanningAlignments.push_back(alignmentResult2);
  BOOST_REQUIRE(!selectJumpContigRNA(candidateAssemblyData, scores));

  // Case-2 is designed here. Out of three jump alignments, alignmentResult3 has max alignment score.
  // And also alignmentResult3 is not Low Quality Jump Alignment as mentioned in
  // test_isLowQualityJumpAlignment. The reason behind high quality is all the jump alignments' score are 1
  // (>0.75).
  JumpAlignmentResult<int> alignmentResult3;
  std::string              testCigar5("35=");
  cigar_to_apath(testCigar5.c_str(), alignmentResult3.align1.apath);  // alignment score is 35
  std::string testCigar6("35=");
  cigar_to_apath(testCigar6.c_str(), alignmentResult3.align2.apath);  // alignment score is 35
  alignmentResult3.score = 70;
  candidateAssemblyData.contigs.resize(3);
  candidateAssemblyData.spanningAlignments.push_back(alignmentResult3);
  BOOST_REQUIRE(selectJumpContigRNA(candidateAssemblyData, scores));
  BOOST_REQUIRE_EQUAL(candidateAssemblyData.bestAlignmentIndex, 2);

  // Case-3 is designed here. Here both alignmentResult3 and alignmentResult4 are
  // supporting high quality spanning SV. Although the score of alignmentResult3 (70)
  // is more than the score of alignmentResult4(64), then also best alignment index
  // is the index of alignmentResult4 as number of supporting reads of contig in
  // alignmentResult4(4) is more than number of supporting reads of contig in
  // alignmentResult3(3).
  JumpAlignmentResult<int> alignmentResult4;
  std::string              testCigar7("32=");
  cigar_to_apath(testCigar7.c_str(), alignmentResult4.align1.apath);  // alignment score is 32
  std::string testCigar8("32=");
  cigar_to_apath(testCigar8.c_str(), alignmentResult4.align2.apath);  // alignment score is 32
  alignmentResult4.score = 64;
  candidateAssemblyData.contigs.resize(4);
  // Contig of alignmentResult1, alignmentResult2 and alignmentResult3 are supported
  // by 3 reads, and contig of alignmentResult4 is supported by 4 reads.
  candidateAssemblyData.contigs[0].supportReads = {1, 2, 3};
  candidateAssemblyData.contigs[1].supportReads = {4, 5, 6};
  candidateAssemblyData.contigs[2].supportReads = {7, 8, 9};
  candidateAssemblyData.contigs[3].supportReads = {10, 11, 12, 13};
  candidateAssemblyData.spanningAlignments.push_back(alignmentResult4);
  BOOST_REQUIRE(selectJumpContigRNA(candidateAssemblyData, scores));
  BOOST_REQUIRE_EQUAL(candidateAssemblyData.bestAlignmentIndex, 3);
}

// Convert jump alignment results into an SVCandidate and add all
// extra data required for VCF output
// Test the following cases:
// 1. Breakend insertion sequence, together with its size, are added to candidate SV.
// 2. If a user wants to output contig sequence in VCF, it will add contig sequence
//    to candidate sv.
// 3. Creation of intervals of candidate SV are explained in test_adjustAssembledBreakend
// 4. Breakend insertion sequence, together with its size, are added to candidate SV. Also
//    breakend deletion size is added to candidate SV.
BOOST_AUTO_TEST_CASE(test_generateRefinedVCFSVCandidateFromJumpAlignment)
{
  // Match score = 1, mismatch penalty = -4,
  // gap opening penalty = -6, gap extension penalty = -2,
  // clipping penalty = -5.
  AlignmentScores<int>     scores(1, -4, -6, -2, -5);
  JumpAlignmentResult<int> alignmentResult1;
  std::string              testCigar1("35=5I30=6D10=3I");             // match, insertion and deletion
  cigar_to_apath(testCigar1.c_str(), alignmentResult1.align1.apath);  // alignment score is 29
  std::string testCigar2("35=");
  cigar_to_apath(testCigar2.c_str(), alignmentResult1.align2.apath);  // alignment score is 35
  alignmentResult1.score = 64;

  JumpAlignmentResult<int> alignmentResult2;
  std::string              testCigar3("35=5X30=6D10=3I");
  cigar_to_apath(testCigar3.c_str(), alignmentResult2.align1.apath);  // alignment score is 25
  std::string testCigar4("35=");
  cigar_to_apath(testCigar4.c_str(), alignmentResult2.align2.apath);  // alignment score is 35
  alignmentResult1.score = 60;

  // Preparing candidate assembly data.
  SVCandidateAssemblyData candidateAssemblyData;
  candidateAssemblyData.spanningAlignments.push_back(alignmentResult1);
  candidateAssemblyData.spanningAlignments.push_back(alignmentResult2);

  // Out of three jump alignments, alignmentResult3 has max alignment score.
  JumpAlignmentResult<int> alignmentResult3;
  alignmentResult3.jumpInsertSize = 5;
  std::string testCigar5("35=");
  cigar_to_apath(testCigar5.c_str(), alignmentResult3.align1.apath);  // alignment score is 35
  std::string testCigar6("35=");
  cigar_to_apath(testCigar6.c_str(), alignmentResult3.align2.apath);  // alignment score is 35
  alignmentResult3.score           = 70;
  alignmentResult3.align1.beginPos = 9;
  alignmentResult3.align2.beginPos = 40;
  alignmentResult3.jumpRange       = 2;
  candidateAssemblyData.spanningAlignments.push_back(alignmentResult3);

  // Preparing dummy contig sequences
  AssembledContig contig1;
  contig1.seq = "TCTATCACCCATCGTACCACTCACGGGAGCTCGATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGT";
  AssembledContig contig2;
  contig2.seq = "AGCTAGTCAGATCGTACCACTCACGGGAGCTCGATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGT";
  AssembledContig contig3;
  contig3.seq = "TGCATGACGTATCGTACCACTCACGGGAGCTCGATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGT";
  candidateAssemblyData.contigs = {contig1, contig2, contig3};
  candidateAssemblyData.bp1ref.seq() =
      "GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGT"
      "ATTTTCGTCTGGGGGGTGTGCACGCGATAGCATTGCGAGACGCTGGA";
  candidateAssemblyData.bp2ref.seq() =
      "TCCAGCGTCTCGCAATGCTATCGCGTGCACACCCCCCAGACGAAAATACCAAATGCATGGAGA"
      "GCTCCCGTGAGTGGTTAATAGGGTGATAGACCTGTGATC";
  // BP2 is stated after 40 bases from the beginning of the reference
  candidateAssemblyData.bp2ref.set_offset(40);
  // Best alignment index as explained in test_selectJumpContigDNA.
  candidateAssemblyData.bestAlignmentIndex = 2;
  SVCandidate candidate;
  candidate.bp1.state = SVBreakendState::RIGHT_OPEN;
  candidate.bp2.state = SVBreakendState::LEFT_OPEN;
  candidate.setPrecise();
  GSCOptions options;
  options.isOutputContig = true;  // output contig sequence to VCF.
  generateRefinedVCFSVCandidateFromJumpAlignment(candidateAssemblyData, candidate, options);
  // After refinement, SV breakpoints are BP1([43, 46)) and BP2([80, 83)).
  // deletion size = 80 - 43 - 1 = 36
  // Insertion size = jump insert size of best alignment = 5.
  std::string       expectedAlignment("5I36D");
  ALIGNPATH::path_t expectedPath;
  cigar_to_apath(expectedAlignment.c_str(), expectedPath);

  // Case-1 is designed where breakend insert sequence size is jump insert size.
  BOOST_REQUIRE_EQUAL(candidate.insertSeq.size(), alignmentResult3.jumpInsertSize);
  // Case-2 is designed. As best jump aligment is alignmentResult3, so candidate contig sequence
  // is the contig sequence of alignmentResult3.
  BOOST_REQUIRE_EQUAL(candidate.contigSeq, contig3.seq);
  // Case-3 is designed here.
  // homologous jump range across the breakend = 2
  // As state is RIGHT_OPEN, range start is alignment end - 1 = 44 - 1 = 43
  // and alignment end = range start + jump range + 1 = 46
  BOOST_REQUIRE_EQUAL(candidate.bp1.interval, GenomeInterval(0, 43, 46));
  // homologous jump range across the breakend = 2
  // As state is LEFT_OPEN, range start is alignment start + ref_offset = 40 + 40 = 80
  // and alignment end = alignment start + jump range + 1 = 83
  BOOST_REQUIRE_EQUAL(candidate.bp2.interval, GenomeInterval(0, 80, 83));
  // Case-4 is designed where expected path is 5I36D.
  BOOST_REQUIRE_EQUAL(candidate.insertAlignment, expectedPath);
}

// Given a SV candidate, Compute a possible assembly.
// Test the following cases:
// For Spanning SV:
// 1. When two breakpoints are overlapping (with extra padding 350) each other, api computes small
//    candidate assembly. This case assumes a single-interval local assembly, this is the most
//    common case for small-scale SVs/indels.
// 2. For large SV or interchromosomal SV, api computes spanning candidate assembly. This case assumes
//    two suspected breakends with a direction to each, most common large scale SV case.
// 3. Above two cases are applied for RNA also.
//    For Complex SV:
// 4. When a SV is complex, api computes small candidate assembly as mentioned in case-1.
//
//  Contigs are generated based on the bam record specified at the top of this file. BamRecords are
//  created in such a way that leading or trailing has 5 mismatches (as 4 is the leading or trailing
//  mismatch threshold)
BOOST_AUTO_TEST_CASE(test_getCandidateAssemblyData)
{
  const bam_header_info bamHeader(buildTestBamHeader());
  TestFilenameMaker     fileMaker1;
  GSCOptions            options;
  // Assembly options
  options.refineOpt.spanningAssembleOpt.minWordLength = 15;
  options.refineOpt.spanningAssembleOpt.maxWordLength = 40;
  options.refineOpt.spanningAssembleOpt.wordStepSize  = 10;
  options.refineOpt.smallSVAssembleOpt.minWordLength  = 3;
  options.refineOpt.smallSVAssembleOpt.maxWordLength  = 9;
  options.refineOpt.smallSVAssembleOpt.wordStepSize   = 3;
  options.alignFileOpt.alignmentFilenames             = {bamFileName};
  options.alignFileOpt.isAlignmentTumor               = {false};
  options.referenceFilename                           = getTestReferenceFilename();
  options.edgeRuntimeFilename                         = fileMaker1.getFilename();
  // Creating stats file.
  TestStatsFileMaker statsFile;
  options.statsFilename = statsFile.getFilename();
  AllSampleReadCounts counts;
  counts.setSampleCount(1);
  SampleReadCounts sample1(counts.getSampleCounts(0));
  sample1.input.evidenceCount.anom         = 6;
  sample1.input.evidenceCount.split        = 0;
  sample1.input.evidenceCount.anomAndSplit = 0;
  sample1.input.evidenceCount.total        = 4;
  counts.getSampleCounts(0).merge(sample1);
  auto edgeTrackerPtr(std::make_shared<EdgeRuntimeTracker>(options.edgeRuntimeFilename));

  // Case-1 is designed here. It is a spanning sv candidate where breakpoints
  // are overlapping (with extra padding 350) each other.
  // As the breakends are overlapping each other using this padding, api treats
  // it as small SV alignment. Here this SV is supported by Read-4, Read-5 and Read-6.
  // Let's denote these 3 reads by some index like 0,1,2 respectively. As two breakpoints
  // are overlapping (using 350 padding) each other, so this padded BP is supported by
  // reads with read indices 0,1 and 2.
  SVCandidateAssemblyRefiner refiner1(options, bamHeader, counts, edgeTrackerPtr);
  SVCandidate                candidate1;
  candidate1.bp1.state    = SVBreakendState::RIGHT_OPEN;
  candidate1.bp1.interval = GenomeInterval(0, 310, 320);
  candidate1.bp2.state    = SVBreakendState::LEFT_OPEN;
  candidate1.bp2.interval = GenomeInterval(0, 330, 350);
  SVCandidateAssemblyData candidateAssemblyData1;
  refiner1.getCandidateAssemblyData(candidate1, false, candidateAssemblyData1);
  BOOST_REQUIRE_EQUAL(candidateAssemblyData1.smallSVAlignments.size(), 1);
  BOOST_REQUIRE_EQUAL(candidateAssemblyData1.contigs.size(), 1);
  // Read-4, Read-5 and Read-6 are supporting the contig.
  BOOST_REQUIRE_EQUAL(candidateAssemblyData1.contigs[0].supportReads.size(), 3);
  // Check read indices based on the current set. As the current set contains only
  // Read-4, Read-5 and Read-6, so their read indices are 0,1 and 2 respectively.
  BOOST_REQUIRE_EQUAL(*(candidateAssemblyData1.contigs[0].supportReads.begin()), 0);
  BOOST_REQUIRE_EQUAL(*(std::next(candidateAssemblyData1.contigs[0].supportReads.begin(), 1)), 1);
  BOOST_REQUIRE_EQUAL(*(std::next(candidateAssemblyData1.contigs[0].supportReads.begin(), 2)), 2);
  BOOST_REQUIRE_EQUAL(candidateAssemblyData1.contigs[0].seq, "AAACCCCCCCCTCCCCCCGCTTCTGGCCTTTTTTTTT");

  // Case-2 designed here. It is a spanning inter-chromosomal SV which is supported by
  // Read-1, Read-2 and Read-3 pairs. Let's denote these 6 reads by some index. Let's say
  // read indices of mate one reads are 0, 1, 2 respectively and read indices of mate two
  // reads are 3, 4, 5 respectively. Then BP1 is supported by reads with read indices 0,1 and 2
  // and BP2 is supported by reads with read indices 3,4 and 5.
  SVCandidate candidate2;
  candidate2.bp1.state    = SVBreakendState::RIGHT_OPEN;
  candidate2.bp1.interval = GenomeInterval(0, 40, 50);
  candidate2.bp2.state    = SVBreakendState::LEFT_OPEN;
  candidate2.bp2.interval = GenomeInterval(1, 65, 75);
  SVCandidateAssemblyData candidateAssemblyData2;
  refiner1.getCandidateAssemblyData(candidate2, false, candidateAssemblyData2);
  // Two contigs are generated for BP1 and BP2.
  BOOST_REQUIRE_EQUAL(candidateAssemblyData2.contigs.size(), 2);
  // Each contig aligns to the reference. So total number of spanning alignment is 2.
  BOOST_REQUIRE_EQUAL(candidateAssemblyData2.spanningAlignments.size(), 2);
  // Read indices {3,4,5} support BP2 contig.
  BOOST_REQUIRE_EQUAL(candidateAssemblyData2.contigs[0].supportReads.size(), 3);
  BOOST_REQUIRE_EQUAL(*(candidateAssemblyData2.contigs[0].supportReads.begin()), 3);
  BOOST_REQUIRE_EQUAL(*(std::next(candidateAssemblyData2.contigs[0].supportReads.begin(), 1)), 4);
  BOOST_REQUIRE_EQUAL(*(std::next(candidateAssemblyData2.contigs[0].supportReads.begin(), 2)), 5);
  BOOST_REQUIRE_EQUAL(candidateAssemblyData2.contigs[0].seq, "AAAAAAACTCATCAGTTGATGATACGCCCGAGCAGAT");
  // Read indices {0,1,2} support BP1 contig.
  BOOST_REQUIRE_EQUAL(candidateAssemblyData2.contigs[1].supportReads.size(), 3);
  BOOST_REQUIRE_EQUAL(*(candidateAssemblyData2.contigs[1].supportReads.begin()), 0);
  BOOST_REQUIRE_EQUAL(*(std::next(candidateAssemblyData2.contigs[1].supportReads.begin(), 1)), 1);
  BOOST_REQUIRE_EQUAL(*(std::next(candidateAssemblyData2.contigs[1].supportReads.begin(), 2)), 2);
  BOOST_REQUIRE_EQUAL(candidateAssemblyData2.contigs[1].seq, "GTCTATCACCCTATTAACCACTCACGGGAGAAAAAGA");

  // Case-3 is designed here. This is for RNA sample and this is same as the above
  // test case.
  SVCandidate candidate3;
  candidate3.bp1.state    = SVBreakendState::RIGHT_OPEN;
  candidate3.bp1.interval = GenomeInterval(0, 40, 50);
  candidate3.bp2.state    = SVBreakendState::LEFT_OPEN;
  candidate3.bp2.interval = GenomeInterval(1, 65, 75);
  SVCandidateAssemblyData candidateAssemblyData3;
  options.isRNA = true;  // This is a RNA sample.
  SVCandidateAssemblyRefiner refiner2(options, bamHeader, counts, edgeTrackerPtr);
  refiner2.getCandidateAssemblyData(candidate3, false, candidateAssemblyData3);
  BOOST_REQUIRE_EQUAL(candidateAssemblyData3.contigs.size(), 2);
  // Read indices {3,4,5} support BP2 contig.
  BOOST_REQUIRE_EQUAL(candidateAssemblyData3.contigs[0].supportReads.size(), 3);
  BOOST_REQUIRE_EQUAL(*(candidateAssemblyData3.contigs[0].supportReads.begin()), 3);
  BOOST_REQUIRE_EQUAL(*(std::next(candidateAssemblyData3.contigs[0].supportReads.begin(), 1)), 4);
  BOOST_REQUIRE_EQUAL(*(std::next(candidateAssemblyData3.contigs[0].supportReads.begin(), 2)), 5);
  BOOST_REQUIRE_EQUAL(candidateAssemblyData3.contigs[0].seq, "AAAAAAACTCATCAGTTGATGATACGCCCGAGCAGAT");
  // Read indices {0,1,2} support BP1 contig.
  BOOST_REQUIRE_EQUAL(candidateAssemblyData3.contigs[1].supportReads.size(), 3);
  BOOST_REQUIRE_EQUAL(*(candidateAssemblyData3.contigs[1].supportReads.begin()), 0);
  BOOST_REQUIRE_EQUAL(*(std::next(candidateAssemblyData3.contigs[1].supportReads.begin(), 1)), 1);
  BOOST_REQUIRE_EQUAL(*(std::next(candidateAssemblyData3.contigs[1].supportReads.begin(), 2)), 2);
  BOOST_REQUIRE_EQUAL(candidateAssemblyData3.contigs[1].seq, "GTCTATCACCCTATTAACCACTCACGGGAGAAAAAGA");

  // Case-4 is designed here. This is a complex SV. This case assumes a single-interval local assembly,
  // this is the most common case for small-scale SVs/indels. For complex SV, smallSVAlignments are
  // expected.
  SVCandidate candidate4;
  candidate4.bp1.state    = SVBreakendState::COMPLEX;
  candidate4.bp1.interval = GenomeInterval(0, 310, 320);
  candidate4.bp2.state    = SVBreakendState::UNKNOWN;
  candidate4.bp2.interval = GenomeInterval(0, 330, 350);
  SVCandidateAssemblyData    candidateAssemblyData4;
  SVCandidateAssemblyRefiner refiner3(options, bamHeader, counts, edgeTrackerPtr);
  refiner3.getCandidateAssemblyData(candidate4, true, candidateAssemblyData4);
  // As this is a complex SV, contigs are generated based on BP1 (BP2 is unknown).
  BOOST_REQUIRE_EQUAL(candidateAssemblyData4.contigs.size(), 1);
  // Each contig aligns to the reference. So total number of alignment is 1.
  BOOST_REQUIRE_EQUAL(candidateAssemblyData4.smallSVAlignments.size(), 1);
  // Read-4, Read-5 and Read-6 are supporting the contig.
  BOOST_REQUIRE_EQUAL(candidateAssemblyData4.contigs[0].supportReads.size(), 3);
  // Read indices based on the current set. As current set contains only
  // Read-4, Read-5 and Read-6, so their read indices are 0,1 and 2 respectively.
  BOOST_REQUIRE_EQUAL(*(candidateAssemblyData4.contigs[0].supportReads.begin()), 0);
  BOOST_REQUIRE_EQUAL(*(std::next(candidateAssemblyData4.contigs[0].supportReads.begin(), 1)), 1);
  BOOST_REQUIRE_EQUAL(*(std::next(candidateAssemblyData4.contigs[0].supportReads.begin(), 2)), 2);
  BOOST_REQUIRE_EQUAL(candidateAssemblyData4.contigs[0].seq, "AAACCCCCCCCTCCCCCCGCTTCTGGCCTTTTTTTTT");
}

// Search for combinations of left and right-side insertion candidates to find a good insertion pair
// Test the following case:
// If there is a large insertion, api creates a new SV candidate around that insertion region.
BOOST_AUTO_TEST_CASE(test_ProcessLargeInsertion)
{
  GSCOptions               options;
  const GlobalAligner<int> largeInsertCompleteAligner(options.refineOpt.largeInsertEdgeAlignScores);
  SVCandidateAssemblyData  candidateAssemblyData;
  std::vector<SVCandidateAssemblyData::SmallAlignmentResultType>& smallSVAlignments =
      candidateAssemblyData.smallSVAlignments;
  smallSVAlignments.resize(2);
  std::string cigar1("5=60I50=");
  smallSVAlignments[0].align.beginPos = 50;
  cigar_to_apath(cigar1.c_str(), smallSVAlignments[1].align.apath);
  smallSVAlignments[0].score = -23;
  std::string cigar2("10=50I50=");
  smallSVAlignments[1].align.beginPos = 50;
  cigar_to_apath(cigar2.c_str(), smallSVAlignments[0].align.apath);
  smallSVAlignments[1].score = -8;

  std::vector<LargeInsertionInfo>& largeInsertInfo = candidateAssemblyData.largeInsertInfo;
  largeInsertInfo.resize(2);
  largeInsertInfo[0].score            = -23;
  largeInsertInfo[0].isRightCandidate = true;
  largeInsertInfo[0].contigOffset     = 5;
  largeInsertInfo[0].refOffset        = 5;
  largeInsertInfo[1].score            = -8;
  largeInsertInfo[1].isLeftCandidate  = true;
  largeInsertInfo[1].contigOffset     = 10;
  largeInsertInfo[1].refOffset        = 10;

  candidateAssemblyData.bp1ref.seq() =
      "AATATTACCCTTAGTTCTCATTAGTAAATGGATGTGAAAGTCAGCCTGACTGCAGTGCTGA"
      "GGGGGCCCCATTTTTTTACCCCAAATATTACCCTTAGTTCTC";
  Assembly& contigs = candidateAssemblyData.contigs;
  contigs.resize(2);
  // These are the left and right contigs. Both contigs are combined and try to align with bp1ref. After
  // aligning the cigar is 53=258I50=.
  contigs[0].seq =
      "TCCTTCCTTTCCTTTCTAAATCTCTCTCTGTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCGCTCTCTCC"
      "TCCCCTAGTTTATCCTGGACTCATGCTGAGCTCAGCAACCCTTG"
      "AGTGCTGAGGGGGCCCCATTTTTTTACCCCAAATATTACCCTTAGTTCTC";
  contigs[1].seq =
      "AATATTACCCTTAGTTCTCATTAGTAAATGGATGTGAAAGTCAGCCTGACTGCATTCCTTTGCATCGTA"
      "TTCCTCAGTCTCTGCGGAAGGCACTGC";
  SVCandidate     candidate;
  std::set<pos_t> excludedPos;
  processLargeInsertion(
      candidate, 0, 0, largeInsertCompleteAligner, {0, 1}, excludedPos, candidateAssemblyData, options);
  std::string       expectedInsertAlignment = "258I";
  ALIGNPATH::path_t expectedAlignPath;
  cigar_to_apath(expectedInsertAlignment.c_str(), expectedAlignPath);
  BOOST_REQUIRE_EQUAL(candidateAssemblyData.svs.size(), 1);
  BOOST_REQUIRE_EQUAL(candidateAssemblyData.svs[0].insertAlignment, expectedAlignPath);
  // As the actual insertion is located after 53 base pair, so that following breakend intervals are created.
  BOOST_REQUIRE_EQUAL(candidateAssemblyData.svs[0].bp1.interval, GenomeInterval(0, 52, 54));
  BOOST_REQUIRE_EQUAL(candidateAssemblyData.svs[0].bp2.interval, GenomeInterval(0, 53, 55));
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
