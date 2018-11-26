//
// Manta - Structural Variant and Indel Caller
// Copyright (c) 2013-2018 Illumina, Inc.
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
#include "test/testAlignmentDataUtil.hh"
#include "htsapi/SimpleAlignment_bam_util.hh"
// test static function in TU:
#include "applications/GenerateSVCandidates/SVCandidateAssemblyRefiner.cpp"


BOOST_AUTO_TEST_SUITE( test_SVRefiner )

BOOST_AUTO_TEST_CASE( test_GetVariantRange )
{
    known_pos_range2 res;

    const std::string seq1("ABCDDABC");
    const std::string seq2("ABCDDDABC");

    {
        // left shifted case:
        const known_pos_range2 seq1Range(3,3);
        const known_pos_range2 seq2Range(3,4);

        // order reflects a deletion
        res = getVariantRange(seq2,seq2Range,seq1,seq1Range);
        BOOST_REQUIRE_EQUAL(res.begin_pos(), 0);
        BOOST_REQUIRE_EQUAL(res.end_pos(), 2);

        // order reflects an insertion
        res = getVariantRange(seq1,seq1Range,seq2,seq2Range);
        BOOST_REQUIRE_EQUAL(res.begin_pos(), 0);
        BOOST_REQUIRE_EQUAL(res.end_pos(), 2);
    }

    {
        // right shifted case:
        const known_pos_range2 seq1Range(5,5);
        const known_pos_range2 seq2Range(5,6);

        // order reflects a deletion
        res = getVariantRange(seq2,seq2Range,seq1,seq1Range);
        BOOST_REQUIRE_EQUAL(res.begin_pos(), -2);
        BOOST_REQUIRE_EQUAL(res.end_pos(), 0);

        // order reflects an insertion
        res = getVariantRange(seq1,seq1Range,seq2,seq2Range);
        BOOST_REQUIRE_EQUAL(res.begin_pos(), -2);
        BOOST_REQUIRE_EQUAL(res.end_pos(), 0);
    }
}

BOOST_AUTO_TEST_CASE_EXPECTED_FAILURES(test_getInsertTrim, 1)

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
BOOST_AUTO_TEST_CASE( test_getInsertTrim )
{
    bam_record bam_record1;
    buildTestBamRecord(bam_record1, 0, 9, 0, 30, 60, 15, "35M5I15M");
    SimpleAlignment simpleAlignment1(getAlignment(bam_record1));
    // segment range is (1, 1)
    known_pos_range2 range1(getInsertTrim(simpleAlignment1.path, std::pair<unsigned, unsigned>(1, 1)));
    // As segment range (1, 1), read segment start will be after 1st element of cigar(35M) and read segment
    // end will be after 2nd element of the cigar. so here read segment start is 35 and read
    // segment end is 40.
    BOOST_REQUIRE_EQUAL(range1.begin_pos(), 35);
    BOOST_REQUIRE_EQUAL(range1.end_pos(), 40);

    bam_record bam_record2;
    buildTestBamRecord(bam_record2, 0, 9, 0, 30, 60, 15, "35M5D5I15M");
    SimpleAlignment simpleAlignment2(getAlignment(bam_record2));
    // segment range is (1, 1)
    known_pos_range2 range2(getInsertTrim(simpleAlignment2.path, std::pair<unsigned, unsigned>(1, 1)));
    // As segment range (1, 1), read segment start will be after 1st element of cigar(35M) and read segment
    // end will be after 2nd element(5D) of the cigar. so here read segment start is 35 and read
    // segment end is 35 (Deletion is not included in the calculation).
    BOOST_REQUIRE_EQUAL(range2.begin_pos(), 35);
    BOOST_REQUIRE_EQUAL(range2.end_pos(), 35);

    bam_record bam_record3;
    buildTestBamRecord(bam_record3, 0, 9, 0, 30, 60, 15, "35M5D15M");
    SimpleAlignment simpleAlignment3(getAlignment(bam_record3));
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
// 3. Size of largest indel among multiple cigar segments
BOOST_AUTO_TEST_CASE( test_getLargestIndelSize )
{
    std::string testCigar1("35M5I30M");
    ALIGNPATH::path_t path1;
    cigar_to_apath(testCigar1.c_str(), path1);
    std::vector<std::pair<unsigned, unsigned>> segments;
    segments.push_back(std::pair<unsigned, unsigned>(1, 1));
    // Here only segment is (1, 1). As there is only one insertion
    // element in this segment, so size of the largest indel is zero.
    BOOST_REQUIRE_EQUAL(getLargestIndelSize(path1, segments), 5);

    segments.clear();
    segments.push_back(std::pair<unsigned, unsigned>(0,0));
    segments.push_back(std::pair<unsigned, unsigned>(2,2));
    // Here segments are (0, 0) and (2, 2). As there is no indel
    // in these two segments of the cigar, so size of the largest
    // indel is zero.
    BOOST_REQUIRE_EQUAL(getLargestIndelSize(path1, segments), 0);

    segments.clear();
    segments.push_back(std::pair<unsigned, unsigned>(0,0));
    segments.push_back(std::pair<unsigned, unsigned>(1,1));
    // Here segments are (0, 0) and (1, 1). As there is no indel
    // in the (0, 0) segment and one insertion element(5I) in the (1, 1)
    // segment of the cigar, so size of the largest indel is 5.
    BOOST_REQUIRE_EQUAL(getLargestIndelSize(path1, segments), 5);

    std::string testCigar2("35M5I30M6D10M");
    ALIGNPATH::path_t path2;
    cigar_to_apath(testCigar2.c_str(), path2);
    segments.clear();
    segments.push_back(std::pair<unsigned, unsigned>(0,0));
    segments.push_back(std::pair<unsigned, unsigned>(1,1));
    segments.push_back(std::pair<unsigned, unsigned>(4,4));
    // Here segments are (0, 0), (1, 1) and (4, 4). Among these 4 segments
    // (1,1) segment has the largest indel and size is 5.
    BOOST_REQUIRE_EQUAL(getLargestIndelSize(path2, segments), 5);

    segments.clear();
    segments.push_back(std::pair<unsigned, unsigned>(0,0));
    segments.push_back(std::pair<unsigned, unsigned>(1,1));
    segments.push_back(std::pair<unsigned, unsigned>(3,4));
    // Here segments are (0, 0), (1, 1) and (4, 4). Among these 4 segments,
    // (3,4) segment has the largest indel and size is 6.
    BOOST_REQUIRE_EQUAL(getLargestIndelSize(path2, segments), 6);
}

// Test the larger indel location in the cigar based on
// minIndel size as indel theshold.
BOOST_AUTO_TEST_CASE( test_getLargeIndelSegments )
{
    std::string testCigar1("35M5I30M6D10M3I");
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

    // Indel threshold = 3. There is three indels with size greater
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
}

// Test the largest insert segment with minInsertThreshold.
BOOST_AUTO_TEST_CASE( test_getLargestInsertSegment )
{
    std::string testCigar1("35M5I30M6D10M2I");
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

    std::string testCigar2("35M5I30M6D10M6I");
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
}

// Test the number of locations when query sequence aligns gaplessly with
// the target sequence provided a threshold of mismatch rate(fraction of mismatch
// allowed).
// If target size is 100 and query size is 50, so from 0 to 49th location
// it will try to align and calculate the mismatch rate.
BOOST_AUTO_TEST_CASE( test_getQuerySeqMatchCount )
{
    std::string targetSeq = "GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGT"
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
BOOST_AUTO_TEST_CASE( test_isLowQualitySmallSVAlignment )
{
    // Match score = 1, mismatch penalty = -4,
    // gap opening penalty = -6, gap extension penalty = -2,
    // clipping penalty = -5.
    AlignmentScores<int> scores(1, -4, -6, -2, -5);
    std::string testCigar1("35=5I30=6D10=3I");
    ALIGNPATH::path_t path1;
    cigar_to_apath(testCigar1.c_str(), path1);
    // According to penalty, alignment score = 29 and optimal score = 83.
    // So score fraction is 29/83 = ~0.34 (<0.75)
    BOOST_REQUIRE(isLowQualitySmallSVAlignment(100, scores, false, false, path1));

    std::string testCigar2("35=5X30=6D10=3I");
    ALIGNPATH::path_t path2;
    cigar_to_apath(testCigar2.c_str(), path2);
    // According to penalty, alignment score = 25 and optimal score = 83.
    // So score fraction is 25/83 = ~0.30 (<0.75)
    BOOST_REQUIRE(isLowQualitySmallSVAlignment(100, scores, false, false, path2));

    std::string testCigar3("35=");
    ALIGNPATH::path_t path3;
    cigar_to_apath(testCigar3.c_str(), path3);
    // According to penalty, alignment score = 35 and optimal score = 35.
    // So score fraction is 35/35 = 1.0 (>=0.75)
    BOOST_REQUIRE(!isLowQualitySmallSVAlignment(100, scores, false, false, path3));

    std::string testCigar4("25=");
    ALIGNPATH::path_t path4;
    cigar_to_apath(testCigar4.c_str(), path4);
    // Here reference projection length is 25 which is less than minAlignRefSpan.
    BOOST_REQUIRE(isLowQualitySmallSVAlignment(100, scores, false, false, path4));

    std::string testCigar5("25=23S");
    ALIGNPATH::path_t path5;
    cigar_to_apath(testCigar5.c_str(), path5);
    // Here clipped path size 25 which is less than minAlignReadLength.
    BOOST_REQUIRE(isLowQualitySmallSVAlignment(100, scores, false, false, path5));

    std::string testCigar6("30=23S");
    ALIGNPATH::path_t path6;
    cigar_to_apath(testCigar6.c_str(), path6);
    // Here reference projection length and clipped path size satisfy the thresholds.
    // Soft-clip will not be considered in score calculation.
    // According to penalty, alignment score = 30 and optimal score = 30.
    // So score fraction is 30/30 = 1.0 (>=0.75)
    BOOST_REQUIRE(!isLowQualitySmallSVAlignment(100, scores, false, false, path6));

    // Following two cases are for complex candidate
    std::string testCigar7("30=");
    ALIGNPATH::path_t path7;
    cigar_to_apath(testCigar7.c_str(), path4);
    // Here reference projection length is 30 which is less than minAlignRefSpan.
    BOOST_REQUIRE(isLowQualitySmallSVAlignment(100, scores, false, true, path7));

    std::string testCigar8("30=23S");
    ALIGNPATH::path_t path8;
    cigar_to_apath(testCigar8.c_str(), path8);
    // Here clipped path size 30 which is less than minAlignReadLength.
    BOOST_REQUIRE(isLowQualitySmallSVAlignment(100, scores, false, true, path8));
}

// Test the following cases:
// 1. Whether simple cigar string is inserted to SV candidate
// 2. Insertion size = size of the candidate insert sequence
// 3. Deletion size = Difference between the two breakpoint start position - 1
// 4. All the above points are valid for indel type sv, otherwise it will not anything.
// Following two schematic diagrams are vaild for indel type SV
//                   BP1               BP2
// Category-1   >>>>>>> -------------- <<<<<<< (BP1 Right-Open and BP2 Left-Open )
//
//                   BP2               BP1
// Category-2   >>>>>>> -------------- <<<<<<< (BP1 Left-Open and BP2 Right-Open )
BOOST_AUTO_TEST_CASE( test_addCigarToSpanningAlignment )
{
    // Category-1 (See docs)
    SVCandidate candidate1;
    candidate1.bp1.state = SVBreakendState::RIGHT_OPEN;
    candidate1.bp2.state = SVBreakendState::LEFT_OPEN;
    candidate1.bp1.interval = GenomeInterval(0,1000,1200);
    candidate1.bp2.interval = GenomeInterval(0, 2000, 2100);
    candidate1.insertSeq = "AGAGCTAGT";
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
    candidate2.bp1.interval = GenomeInterval(0,2000,2100);
    candidate2.bp2.interval = GenomeInterval(0, 1000, 1200);
    candidate2.bp1.state = SVBreakendState::LEFT_OPEN;
    candidate2.bp2.state = SVBreakendState::RIGHT_OPEN;
    candidate2.insertSeq = "AGAGCTAGT";
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
    candidate3.bp1.state = SVBreakendState::RIGHT_OPEN;
    candidate3.bp2.state = SVBreakendState::RIGHT_OPEN;
    candidate3.bp1.interval = GenomeInterval(0,1000,1200);
    candidate3.bp2.interval = GenomeInterval(0, 2000, 2100);
    candidate3.insertSeq = "AGAGCTAGT";
    addCigarToSpanningAlignment(candidate3);
    BOOST_REQUIRE_EQUAL(candidate3.insertAlignment.size(), 0);
}

// Test whether max score contiguous path(starting from the beginning) of alignment
// satisfies the following criteria:
// 1. Minimum contig length in the max score contiguous path = 40
// 2. Minimum reference projection length in the max score contiguous path = 40
// 3. Fraction of alignment score of max score path should be greater than 0.75 where
//    fraction of alignment score is calculated as
//    alignment_score / (clipped_contig_length * match score)
// 4. Remaining portion means (read length -contig length) should be greater than 40.
// For example:  If a cigar is 35=5I30=6D10=3I, then max score contiguous path
// is 35=5I30=.
BOOST_AUTO_TEST_CASE( test_isLargeInsertSegment )
{
    // Match score = 1, mismatch penalty = -4,
    // gap opening penalty = -6, gap extension penalty = -2,
    // clipping penalty = -5.
    AlignmentScores<int> scores(1, -4, -6, -2, -5);
    AlignerBase<int> alignerBase(scores);

    // contig length in the max score contiguous path
    unsigned contigOffset1 = 0;
    // Reference projection length in the max score contiguous path
    unsigned refoffset1 = 0;
    // alignment score in the max score contiguous path
    int alignmentScore1 = 0;
    std::string testCigar1("35=5I30=6D10=3I");
    ALIGNPATH::path_t path1;
    cigar_to_apath(testCigar1.c_str(), path1);
    // According to 35=5I30=6D10=3I, max score contiguous path is 35=5I30=
    // But here remaining portion  = 83 - 70 = 13 which is less than 40.
    // So the api will return false.
    BOOST_REQUIRE(!isLargeInsertSegment(alignerBase, path1, contigOffset1, refoffset1, alignmentScore1));
    // Length of the offsets according to 35=5I30=
    BOOST_REQUIRE_EQUAL(contigOffset1, 70);
    BOOST_REQUIRE_EQUAL(refoffset1, 65);
    BOOST_REQUIRE_EQUAL(alignmentScore1, 49);

    unsigned contigOffset2;
    unsigned refoffset2;
    int alignmentScore2;
    std::string testCigar2("100=");
    ALIGNPATH::path_t path2;
    cigar_to_apath(testCigar2.c_str(), path2);
    // According to 100=, max score contiguous path is 100=
    // But here remaining portion  = 100 - 100 = 0 which is less than 40.
    // So the api will return false.
    BOOST_REQUIRE(!isLargeInsertSegment(alignerBase, path2, contigOffset2, refoffset2, alignmentScore2));
    // Length of the offsets according to 100=
    BOOST_REQUIRE_EQUAL(contigOffset2, 100);
    BOOST_REQUIRE_EQUAL(refoffset2, 100);
    BOOST_REQUIRE_EQUAL(alignmentScore2, 100);

    unsigned contigOffset3;
    unsigned refoffset3;
    int alignmentScore3;
    std::string testCigar3("200=40X"); // maxpath is 200=
    ALIGNPATH::path_t path3;
    cigar_to_apath(testCigar3.c_str(), path3);
    // According to 200=40X, max score contiguous path is 200=
    // But here remaining portion  = 240 - 200 = 40 which is fine.
    // And the fraction of alignment score = 200 / 200 = 1 (>0.75).
    BOOST_REQUIRE(isLargeInsertSegment(alignerBase, path3, contigOffset3, refoffset3, alignmentScore3));
    // Length of the offsets according to 200=40X
    BOOST_REQUIRE_EQUAL(contigOffset3, 200);
    BOOST_REQUIRE_EQUAL(refoffset3, 200);
    BOOST_REQUIRE_EQUAL(alignmentScore3, 200);

    unsigned contigOffset4;
    unsigned refoffset4;
    int alignmentScore4;
    std::string testCigar4("49X200=40X");
    ALIGNPATH::path_t path4;
    cigar_to_apath(testCigar4.c_str(), path4);
    // According to 49X200=40X, max score contiguous path is 49X200=
    // But here remaining portion  = 289- 249 = 40 which is fine.
    // And the fraction of alignment score = 4 / 289 = ~0.01 (<0.75).
    BOOST_REQUIRE(!isLargeInsertSegment(alignerBase, path4, contigOffset4, refoffset4, alignmentScore4));
    // Length of the offsets according to 200=40X
    BOOST_REQUIRE_EQUAL(contigOffset4, 249);
    BOOST_REQUIRE_EQUAL(refoffset4, 249);
    BOOST_REQUIRE_EQUAL(alignmentScore4, 4);

    std::string testCigar5("25=1I14=");
    ALIGNPATH::path_t path5;
    cigar_to_apath(testCigar5.c_str(), path5);
    // According to 25=1I14=, max score contiguous path is 25=1I14=
    // But here reference projection length is 39 which is less than 40.
    // So the api will return false.
    BOOST_REQUIRE(!isLargeInsertSegment(alignerBase, path5, contigOffset3, refoffset3, alignmentScore3));

    std::string testCigar6("25=1D14=");
    ALIGNPATH::path_t path6;
    cigar_to_apath(testCigar6.c_str(), path6);
    // According to 25=1D14=, max score contiguous path is 25=1D14=
    // But here contig length is 39 which is less than 40.
    // So the api will return false.
    BOOST_REQUIRE(!isLargeInsertSegment(alignerBase, path6, contigOffset3, refoffset3, alignmentScore3));
}

// Test whether an alignment satisfies the test_isLargeInsertSegment
// either from the beginning or from the end or both.
BOOST_AUTO_TEST_CASE( test_isLargeInsertAlignment )
{
    // Match score = 1, mismatch penalty = -4,
    // gap opening penalty = -6, gap extension penalty = -2,
    // clipping penalty = -5.
    AlignmentScores<int> scores(1, -4, -6, -2, -5);
    GlobalAligner <int> globalAligner(scores);
    LargeInsertionInfo largeInsertionInfo1;
    std::string testCigar1("200=40X");
    ALIGNPATH::path_t path1;
    cigar_to_apath(testCigar1.c_str(), path1);
    // From start it is satisfied test_isLargeInsertSegment but from the end
    // it is not satisfied test_isLargeInsertSegment.
    BOOST_REQUIRE(isLargeInsertAlignment(globalAligner, path1, largeInsertionInfo1));
    BOOST_REQUIRE(largeInsertionInfo1.isLeftCandidate);
    BOOST_REQUIRE(!largeInsertionInfo1.isRightCandidate);

    LargeInsertionInfo largeInsertionInfo2;
    std::string testCigar2("40X200=");
    ALIGNPATH::path_t path2;
    cigar_to_apath(testCigar2.c_str(), path2);
    // From start it is not satisfied test_isLargeInsertSegment but from the end
    // it is satisfied test_isLargeInsertSegment.
    BOOST_REQUIRE(isLargeInsertAlignment(globalAligner, path2, largeInsertionInfo2));
    BOOST_REQUIRE(!largeInsertionInfo2.isLeftCandidate);
    BOOST_REQUIRE(largeInsertionInfo2.isRightCandidate);

    LargeInsertionInfo largeInsertionInfo3;
    std::string testCigar3("150="); // maxpath is 200=
    ALIGNPATH::path_t path3;
    cigar_to_apath(testCigar3.c_str(), path3);
    // From start and end it is not satisfied test_isLargeInsertSegment.
    BOOST_REQUIRE(!isLargeInsertAlignment(globalAligner, path3, largeInsertionInfo3));
    BOOST_REQUIRE(!largeInsertionInfo3.isLeftCandidate);
    BOOST_REQUIRE(!largeInsertionInfo3.isRightCandidate);
}

// Test whether an alignment satisfies test_isLargeInsertSegment
// in a partular segment region.
BOOST_AUTO_TEST_CASE( test_isFinishedLargeInsertAlignment )
{
    AlignmentScores<int> scores(1, -4, -6, -1, -5);
    GlobalAligner <int> globalAligner(scores);
    std::string testCigar1("40I350=40I");
    ALIGNPATH::path_t path1;
    cigar_to_apath(testCigar1.c_str(), path1);
    BOOST_REQUIRE(isFinishedLargeInsertAlignment(globalAligner, path1, std::pair<unsigned, unsigned>(0,2), 0));
    BOOST_REQUIRE(!isFinishedLargeInsertAlignment(globalAligner, path1, std::pair<unsigned, unsigned>(0,2), 1));
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
BOOST_AUTO_TEST_CASE( test_translateMaskedPos )
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
    exclusion_block block3(22, 10, 3);
    std::vector<exclusion_block> blocks = {block1, block2, block3};
    // Original coordinate = 5 + (10-3) + (10-3) = 19 (Based on
    // the above diagram)
    BOOST_REQUIRE_EQUAL(translateMaskedPos(blocks, 5), 19);
}

// Test the simple reference coordinates that can be reported in the output vcf
// Test the following cases:
// 1. In case of sv breakend state is RIGHT_OPEN (coordinates should be around alignment end)
// 2. In case of sv breaked state is LEFT_OPEN (coordinate should be around in alignment start)
BOOST_AUTO_TEST_CASE( test_adjustAssembledBreakend )
{
    reference_contig_segment referenceContigSegment;
    referenceContigSegment.seq() = "GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGT"
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
// 2. Fraction of alignment score should be greater than 0.75 where
//    fraction of alignment score is calculated as
//    alignment_score / (clipped_read_length * match score)
// 3. Soft clip should not be considered.
BOOST_AUTO_TEST_CASE( test_isLowQualitySpanningSVAlignment )
{
    // Match score = 1, mismatch penalty = -4,
    // gap opening penalty = -6, gap extension penalty = -2,
    // clipping penalty = -5.
    AlignmentScores<int> scores(1, -4, -6, -2, -5);
    std::string testCigar1("35=5I30=6D10=3I"); // match, insertion and deletion
    ALIGNPATH::path_t path1;
    cigar_to_apath(testCigar1.c_str(), path1);
    // According to penalty, alignment score = 29 and optimal score = 83.
    // So score fraction is 29/83 = ~0.34 (<0.75)
    BOOST_REQUIRE(isLowQualitySpanningSVAlignment(100, scores, false, false, path1));

    std::string testCigar2("35=5X30=6D10=3I");// match, mismatch, insertion and deletion
    ALIGNPATH::path_t path2;
    cigar_to_apath(testCigar2.c_str(), path2);
    // According to penalty, alignment score = 29 and optimal score = 83.
    // So score fraction is 25/83 = ~0.30 (<0.75)
    BOOST_REQUIRE(isLowQualitySpanningSVAlignment(100, scores, false, false, path2));

    std::string testCigar3("35=");
    ALIGNPATH::path_t path3;
    cigar_to_apath(testCigar3.c_str(), path3);
    // According to penalty, alignment score = 35 and optimal score = 35.
    // So score fraction is 35/35 = 1.0 (>=0.75)
    BOOST_REQUIRE(!isLowQualitySpanningSVAlignment(100, scores, false, false, path3));

    std::string testCigar4("25=");
    ALIGNPATH::path_t path4;
    cigar_to_apath(testCigar4.c_str(), path4);
    // Here read size 25 which is less than minAlignReadLength.
    BOOST_REQUIRE(isLowQualitySpanningSVAlignment(100, scores, false, false, path4));

    std::string testCigar5("25=23S");
    ALIGNPATH::path_t path5;
    cigar_to_apath(testCigar5.c_str(), path5);
    // Here clipped path size 25 which is less than minAlignReadLength.
    BOOST_REQUIRE(isLowQualitySpanningSVAlignment(100, scores, false, false, path5));

    std::string testCigar6("30=23S");
    ALIGNPATH::path_t path6;
    cigar_to_apath(testCigar6.c_str(), path6);
    // According to penalty, non-clipped alignment score = 30 and optimal score = 30.
    // So score fraction is 30/30 = 1.0 (>=0.75)
    BOOST_REQUIRE(!isLowQualitySpanningSVAlignment(100, scores, false, false, path6));

    // Following two cases are for RNA
    std::string testCigar7("15=");
    ALIGNPATH::path_t path7;
    cigar_to_apath(testCigar7.c_str(), path4);
    // Here read size 25 which is less than minAlignReadLength.
    BOOST_REQUIRE(isLowQualitySpanningSVAlignment(100, scores, false, true, path7));

    std::string testCigar8("15=23S");
    ALIGNPATH::path_t path8;
    cigar_to_apath(testCigar8.c_str(), path8);
    // Here clipped path size 15 which is less than minAlignReadLength.
    BOOST_REQUIRE(isLowQualitySpanningSVAlignment(100, scores, false, true, path8));
}

// Test test_adjustAssembledBreakend when jump alignment results to sv candidates
BOOST_AUTO_TEST_CASE( test_generateRefinedSVCandidateFromJumpAlignment )
{
    // Temporary assembly data
    SVCandidateAssemblyData svCandidateAssemblyData;
    // index of the alignment
    svCandidateAssemblyData.bestAlignmentIndex = 0;
    // Created two alignments with jump size = 2
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
    jumpAlignmentResultType.jumpRange = 2;
    svCandidateAssemblyData.spanningAlignments.push_back(jumpAlignmentResultType);
    SVCandidate svCandidate1;
    svCandidate1.bp1.state = SVBreakendState::RIGHT_OPEN;
    svCandidate1.bp2.state = SVBreakendState::RIGHT_OPEN;
    generateRefinedSVCandidateFromJumpAlignment(svCandidateAssemblyData, svCandidate1);
    // results should be according to test_adjustAssembledBreakend
    BOOST_REQUIRE_EQUAL(svCandidate1.bp1.interval.range.begin_pos(), 499);
    BOOST_REQUIRE_EQUAL(svCandidate1.bp1.interval.range.end_pos(), 502);
    BOOST_REQUIRE_EQUAL(svCandidate1.bp2.interval.range.begin_pos(), 620);
    BOOST_REQUIRE_EQUAL(svCandidate1.bp2.interval.range.end_pos(), 623);

    // Created two alignments with jump size = 0
    alignment1.beginPos = 397;
    std::string testCigar3("92=");
    cigar_to_apath(testCigar3.c_str(), alignment1.apath);
    alignment2.beginPos = 510;
    std::string testCigar4("110=75I1D2=");
    cigar_to_apath(testCigar4.c_str(), alignment2.apath);
    jumpAlignmentResultType.align1 = alignment1;
    jumpAlignmentResultType.align2 = alignment2;
    jumpAlignmentResultType.jumpRange = 0;
    svCandidateAssemblyData.spanningAlignments.clear();
    svCandidateAssemblyData.spanningAlignments.push_back(jumpAlignmentResultType);
    SVCandidate svCandidate2;
    svCandidate2.bp1.state = SVBreakendState::RIGHT_OPEN;
    svCandidate2.bp2.state = SVBreakendState::LEFT_OPEN;
    generateRefinedSVCandidateFromJumpAlignment(svCandidateAssemblyData, svCandidate2);
    // results should be according to test_adjustAssembledBreakend
    BOOST_REQUIRE_EQUAL(svCandidate2.bp1.interval.range.begin_pos(), 488);
    BOOST_REQUIRE_EQUAL(svCandidate2.bp1.interval.range.end_pos(), 489);
    BOOST_REQUIRE_EQUAL(svCandidate2.bp2.interval.range.begin_pos(), 510);
    BOOST_REQUIRE_EQUAL(svCandidate2.bp2.interval.range.end_pos(), 511);
}

// Test whether jump alignment QC fails or not.
// Test the following cases:
// 1. Cigar should be non-empty
// 2. Minimum read length should be 20
// 3. QC should fail if one of the alignment fails one of the above two condidtions.
BOOST_AUTO_TEST_CASE( test_isJumpAlignmentQCFail )
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

BOOST_AUTO_TEST_SUITE_END()
