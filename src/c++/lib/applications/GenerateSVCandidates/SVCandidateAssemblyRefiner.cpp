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

/// \file
/// \author Ole Schulz-Trieglaff
/// \author Chris Saunders
/// \author Xiaoyu Chen
/// \author Naoki Nariai
///

#include "SVCandidateAssemblyRefiner.hpp"

#include <iostream>
#include <unordered_set>

#include "alignment/AlignmentScoringUtil.hpp"
#include "alignment/AlignmentUtil.hpp"
#include "blt_util/align_path.hpp"
#include "blt_util/log.hpp"
#include "blt_util/seq_printer.hpp"
#include "blt_util/seq_util.hpp"
#include "common/Exceptions.hpp"
#include "htsapi/samtools_fasta_util.hpp"
#include "manta/SVCandidateUtil.hpp"
#include "manta/SVReferenceUtil.hpp"

//#define DEBUG_REFINER
//#define DEBUG_CONTIG
//#define DEBUG_KMER
//#define DEBUG_VARR

#ifdef DEBUG_REFINER
#include "blt_util/seq_printer.hpp"
#endif

/// process assembly/align info into simple reference coordinates that can be reported in the output vcf:
///
/// \param[in] isAlign1 if true, this breakend was aligned first by the jump aligner, and therefore
/// left-aligned (if fwd) or right-aligned (if rev)
///
/// \param[in] jumpRange homologous range across the breakend
///
static void adjustAssembledBreakend(
    const Alignment&                align,
    const bool                      isAlign1,
    const unsigned                  jumpRange,
    const reference_contig_segment& ref,
    const bool                      isReversed,
    SVBreakend&                     bp)
{
  const pos_t bpBeginOffset(getAlignBeginOffset(align, ref.seq().size(), isReversed));
  const pos_t bpEndOffset(getAlignEndOffset(align, ref.seq().size(), isReversed));

  const bool isBpAtAlignEnd(bp.state == SVBreakendState::RIGHT_OPEN);

  const pos_t bpBreakendOffset(isBpAtAlignEnd ? (bpEndOffset - 1) : bpBeginOffset);
  const pos_t bpBreakendPos(ref.get_offset() + bpBreakendOffset);

  const bool isLeftAligned(isAlign1 == isBpAtAlignEnd);

  known_pos_range2& range(bp.interval.range);

  if (isLeftAligned) {
    range.set_begin_pos(bpBreakendPos);
    range.set_end_pos(bpBreakendPos + static_cast<pos_t>(jumpRange) + 1);
  } else {
    range.set_begin_pos(bpBreakendPos - static_cast<pos_t>(jumpRange));
    range.set_end_pos(bpBreakendPos + 1);
  }
}

/// \brief Determine if a candidate spanning SV alignment should be filtered due to low quality
///
/// \param[in] maxQCRefSpan Longest flanking sequence length considered for the high quality requirement
/// \param[in] alignmentScores Alignment scores to use for quality assessment of \p input_apath
/// \return True if the alignment path is low quality and should be filtered out
static bool isLowQualitySpanningSVAlignment(
    const unsigned              maxQCRefSpan,
    const AlignmentScores<int>& alignmentScores,
    const bool                  isLeadingPath,
    const bool                  isRNA,
    const ALIGNPATH::path_t&    input_apath)
{
  /// Require min length of each contig sub-alignment even after off-reference clipping
  unsigned minAlignReadLength(30);

  /// Require min fraction of optimal score in each contig sub-alignment
  static const float minScoreFrac(0.75);
  if (isRNA) {
    minAlignReadLength = 20;
  }

  ALIGNPATH::path_t apath(input_apath);

#ifdef DEBUG_REFINER
  log_os << __FUNCTION__ << ": apath: " << apath << " ; maxRefSpan: " << maxQCRefSpan << "\n";
#endif

  // prepare apath by orienting it always going forward from the breakend and limiting the length to
  // the first maxQCRefSpan ref bases covered:
  //
  if (isLeadingPath) {
    std::reverse(apath.begin(), apath.end());
  }

  apath_limit_ref_length(maxQCRefSpan, apath);

  const unsigned readSize(apath_read_length(apath));
  const unsigned clipSize(apath_soft_clip_right_size(apath));

  assert(clipSize <= readSize);

  const unsigned clippedReadSize(readSize - clipSize);

  if (clippedReadSize < minAlignReadLength) {
#ifdef DEBUG_REFINER
    log_os
        << __FUNCTION__
        << ": Rejecting highest scoring contig sub-alignment. Sub-alignment read length after clipping is: "
        << clippedReadSize << " min size is: " << minAlignReadLength << "\n";
#endif
    return true;
  }

  const int nonClipScore(std::max(0, getPathScore(alignmentScores, apath)));
#ifdef DEBUG_REFINER
  log_os << __FUNCTION__ << ": clippedReadSize=" << clippedReadSize << " matchScore=" << alignmentScores.match
         << " nonClipScore=" << nonClipScore << "\n";
#endif

  const int optimalScore(clippedReadSize * alignmentScores.match);
  assert(optimalScore > 0);

  const float scoreFrac(static_cast<float>(nonClipScore) / static_cast<float>(optimalScore));

  const bool isLowQualityAlignmentScore(scoreFrac < minScoreFrac);

#ifdef DEBUG_REFINER
  if (isLowQualityAlignmentScore) {
    log_os << __FUNCTION__
           << ": Rejecting highest scoring contig sub-alignment. Fraction of optimal alignment score is: "
           << scoreFrac << " minScoreFrac: " << minScoreFrac << "\n";
  }
#endif

  return isLowQualityAlignmentScore;
}

/// Identify indels larger than \p minSize in the input alignment path (\p apath)
///
/// \param[in] apath Input alignment which will be searched for large indels
///
/// \param[in] minSize Minimum qualifying indel size
///
/// \param[out] segments Return the indices in apath for each qualifying indel. Each qualifying indel has a
/// (start,end) index pair to account for combined insetion/deletion (ie. 'complex') events.
static void getLargeIndelSegments(
    const ALIGNPATH::path_t&                    apath,
    const unsigned                              minSize,
    std::vector<std::pair<unsigned, unsigned>>& segments)
{
  using namespace ALIGNPATH;

  bool     isInSegment(false);
  bool     isCandidate(false);
  unsigned segmentStart(0);

  segments.clear();

  const unsigned as(apath.size());
  for (unsigned i(0); i < as; ++i) {
    const path_segment& ps(apath[i]);

    if ((ps.type == DELETE) || (ps.type == INSERT)) {
      if (ps.length >= minSize) isCandidate = true;
      if (!isInSegment) segmentStart = i;
      isInSegment = true;
    } else {
      if (isCandidate) {
        assert(i > 0);
        segments.push_back(std::make_pair(segmentStart, (i - 1)));
      }
      isInSegment = false;
      isCandidate = false;
    }
  }

  if (isCandidate) {
    assert(as > 0);
    segments.push_back(std::make_pair(segmentStart, (as - 1)));
  }
}

static unsigned getLargestIndelSize(
    const ALIGNPATH::path_t& apath, const std::vector<std::pair<unsigned, unsigned>>& segments)
{
  unsigned largestSize(0);

  typedef std::pair<unsigned, unsigned> segment_t;
  for (const segment_t& seg : segments) {
    for (unsigned i(seg.first); i <= seg.second; i++) {
      const ALIGNPATH::path_segment& ps(apath[i]);
      if ((ps.type == ALIGNPATH::DELETE) || (ps.type == ALIGNPATH::INSERT)) {
        if (ps.length > largestSize) largestSize = ps.length;
      }
    }
  }

  return largestSize;
}

/// identify the single largest insert segment, if one exists above minSize:
///
static void getLargestInsertSegment(
    const ALIGNPATH::path_t&                    apath,
    const unsigned                              minSize,
    std::vector<std::pair<unsigned, unsigned>>& segments)
{
  using namespace ALIGNPATH;

  bool     isInSegment(false);
  bool     isCandidate(false);
  unsigned segmentStart(0);

  bool                          isMaxSegment(false);
  unsigned                      maxSegmentSize(minSize);
  std::pair<unsigned, unsigned> maxSegment;

  segments.clear();

  const unsigned as(apath.size());
  for (unsigned i(0); i < as; ++i) {
    const path_segment& ps(apath[i]);

    if ((ps.type == DELETE) || (ps.type == INSERT)) {
      if (ps.type == INSERT) {
        if (ps.length >= maxSegmentSize) {
          isMaxSegment   = true;
          maxSegmentSize = ps.length;
          isCandidate    = true;
        }
      }
      if (!isInSegment) segmentStart = i;
      isInSegment = true;
    } else {
      if (isCandidate) {
        assert(i > 0);
        maxSegment = std::make_pair(segmentStart, (i - 1));
      }
      isInSegment = false;
      isCandidate = false;
    }
  }

  if (isCandidate) {
    assert(as > 0);
    maxSegment = std::make_pair(segmentStart, (as - 1));
  }

  if (isMaxSegment) {
    segments.push_back(maxSegment);
  }
}

/// add simple cigar string to spanning alignments for the subset of cases (insertions and deletions) where
/// this is possible
///
/// note that we may not always print this out, even though we compute the cigar here -- this is dependent on
/// the output file format and conventions related to variant size, precision, etc.
///
static void addCigarToSpanningAlignment(SVCandidate& sv)
{
  const SV_TYPE::index_t svType(getSVType(sv));

  if (svType != SV_TYPE::INDEL) return;

  const bool isBp1First(sv.bp1.interval.range.begin_pos() <= sv.bp2.interval.range.begin_pos());

  const SVBreakend& bpA(isBp1First ? sv.bp1 : sv.bp2);
  const SVBreakend& bpB(isBp1First ? sv.bp2 : sv.bp1);

  const pos_t deleteSize((bpB.interval.range.begin_pos() - bpA.interval.range.begin_pos()) - 1);
  const pos_t insertSize(sv.insertSeq.size());

  assert(deleteSize >= 0);
  assert(insertSize || deleteSize);

  // follow convention from non-spanning alignments of xIyD in case of complex variant:
  if (insertSize) {
    sv.insertAlignment.emplace_back(ALIGNPATH::INSERT, insertSize);
  }
  if (deleteSize) {
    sv.insertAlignment.emplace_back(ALIGNPATH::DELETE, deleteSize);
  }
}

/// \brief Determine if a candidate small SV alignment should be filtered due to low quality
///
/// \param[in] maxQCRefSpan Longest flanking sequence length considered for the high quality requirement
/// \param[in] alignmentScores Alignment scores to use for quality assessment of \p apath
/// \return True if the alignment path is low quality and should be filtered out
static bool isLowQualitySmallSVAlignment(
    const unsigned              maxQCRefSpan,
    const AlignmentScores<int>& alignmentScores,
    const bool                  isLeadingPath,
    const bool                  isComplex,
    ALIGNPATH::path_t&          apath)
{
  // Min reference length for alignment:
  static const unsigned minAlignRefSpanSimple(30);

  // Min length of alignment after off-reference clipping:
  static const unsigned minAlignReadLengthSimple(30);

  // Min reference length for alignment:
  static const unsigned minAlignRefSpanComplex(35);

  // Min length of alignment after off-reference clipping:
  static const unsigned minAlignReadLengthComplex(35);

  // Min reference length for alignment
  static const float minScoreFrac(0.75);

  const unsigned minAlignRefSpan(isComplex ? minAlignRefSpanComplex : minAlignRefSpanSimple);
  const unsigned minAlignReadLength(isComplex ? minAlignReadLengthComplex : minAlignReadLengthSimple);

  // prepare apath by orienting it always going forward from the breakend and limiting the length to
  // the first maxQCRefSpan ref bases covered:
  //
  if (isLeadingPath) {
    std::reverse(apath.begin(), apath.end());
  }

  apath_limit_ref_length(maxQCRefSpan, apath);

  const unsigned refSize(apath_ref_length(apath));
  if (refSize < minAlignRefSpan) {
    return true;
  }

  const unsigned pathSize(apath_read_length(apath));
  const unsigned clipSize(apath_soft_clip_right_size(apath));

  assert(clipSize <= pathSize);

  const unsigned clippedPathSize(pathSize - clipSize);

  if (clippedPathSize < minAlignReadLength) {
#ifdef DEBUG_REFINER
    log_os << "Rejecting highest scoring contig sub-alignment. isFirst?: " << isFirstRead
           << ". Sub-alignmnet read length after clipping is: " << clippedPathSize
           << " min size is: " << minAlignReadLength << "\n";
#endif
    return true;
  }

  const int nonClipScore(std::max(0, getPathScore(alignmentScores, apath)));
  const int optimalScore(clippedPathSize * alignmentScores.match);

  const float scoreFrac(static_cast<float>(nonClipScore) / static_cast<float>(optimalScore));
  const bool  isLowQualityAlignmentScore(scoreFrac < minScoreFrac);

#ifdef DEBUG_REFINER
  if (isLowQualityAlignmentScore) {
    log_os << "Rejecting highest scoring contig sub-alignment. isFirst?: " << isFirstRead
           << ". Fraction of optimal alignment score is: " << scoreFrac << " minScoreFrac: " << minScoreFrac
           << "\n";
  }
#endif

  return isLowQualityAlignmentScore;
}

/// Return the number of locations where the \p querySeq sequence matches the \p targetSeq sequence, allowing
/// for a maximum mismatch rate, \p mismatchRate, for each sequence match.
///
static int getQuerySeqMatchCount(
    const std::string& targetSeq, const std::string& querySeq, const float maxMismatchRate)
{
  unsigned       querySeqMatchCount(0);
  const unsigned querySize(querySeq.size());
  const unsigned targetSize(targetSeq.size());

  if (querySize > targetSize) return querySeqMatchCount;

  // set the scanning start & end to make sure the candidate window is overlapping the breakpoint
  const unsigned scanStart = 0;
  const unsigned scanEnd   = targetSize - querySize;

  for (unsigned i(scanStart); i <= scanEnd; i++) {
    unsigned mismatches = 0;
    for (unsigned j(0); j < querySize; j++) {
      if ((querySeq[j] != targetSeq[i + j]) || (querySeq[j] == 'N')) mismatches++;
    }

    if (float(mismatches) / float(querySize) <= maxMismatchRate) {
      querySeqMatchCount++;
    }
  }

  return querySeqMatchCount;
}

/// Find whether the input single-node contig alignment contains an interesting variant above the minimum size
/// and
/// (2) passes QC otherwise (appropriate flanking regions, etc)
///
/// \param[in] maxQCRefSpan Longest flanking sequence length considered for the high quality qc requirement
/// \param[out] candidateSegments Segment indices of all qualifying SV candidates are returned in this
/// structure
///
/// \return True if a qualifying candidate variant is found in the input alignment
///
static bool findCandidateVariantsFromComplexSVContigAlignment(
    const unsigned                              maxQCRefSpan,
    const AlignmentScores<int>&                 alignmentScores,
    const Alignment&                            align,
    const std::string&                          contigSeq,
    const std::string&                          refSeq,
    const unsigned                              minCandidateIndelSize,
    std::vector<std::pair<unsigned, unsigned>>& candidateSegments)
{
  using namespace ALIGNPATH;

  const path_t& apath(align.apath);

  // (1) identify all indels above minimum size:
  //
  getLargeIndelSegments(apath, minCandidateIndelSize, candidateSegments);

  // escape if there are no indels above the minimum size
  if (candidateSegments.empty()) return false;

  // Here, 'complex' means any indel besides a single simple insertion or deletion:
  const bool isComplex(
      (candidateSegments.size() > 1) || (candidateSegments[0].first != candidateSegments[0].second));

  // loop through possible leading segments until a clean one is found:
  //
  while (true) {
    // test quality of alignment segments surrounding the variant region:
    const unsigned firstCandIndelSegment(candidateSegments.front().first);
    path_t         leadingPath(apath.begin(), apath.begin() + firstCandIndelSegment);

    static const bool isLeadingPath(true);
    if (!isLowQualitySmallSVAlignment(maxQCRefSpan, alignmentScores, isLeadingPath, isComplex, leadingPath)) {
      break;
    }

    // escape if this was the last segment
    if (1 == candidateSegments.size()) return false;

    candidateSegments.erase(candidateSegments.begin());
  }

  // loop through possible trailing segments until a clean one is found:
  //
  while (true) {
    // test quality of alignment segments surrounding the variant region:
    const unsigned lastCandIndelSegment(candidateSegments.back().second);
    path_t         trailingPath(apath.begin() + lastCandIndelSegment + 1, apath.end());

    static const bool isLeadingPath(false);
    if (!isLowQualitySmallSVAlignment(
            maxQCRefSpan, alignmentScores, isLeadingPath, isComplex, trailingPath)) {
      break;
    }

    // escape if this was the last segment
    if (1 == candidateSegments.size()) return false;

    candidateSegments.pop_back();
  }

  // Filter out candidates with ambiguous contig alignments
  //
  // For the contig segment on each side of the candidate SV, search for more than one position on the
  // reference which the contig segment can be mapped to with 95% identity or higher. If found, then filter
  // out this candidate due to ambiguous contig alignment.
  //
  {
    // TODO: iterate on all segments
    const path_t apathTillSvStart(apath.begin(), apath.begin() + candidateSegments.front().first);
    const path_t apathTillSvEnd(apath.begin(), apath.begin() + candidateSegments.back().second + 1);

    const int         leftSize    = apath_read_length(apathTillSvStart);
    const int         endPos      = apath_read_length(apathTillSvEnd);
    const int         rightSize   = contigSeq.length() - endPos;
    const std::string leftContig  = contigSeq.substr(0, leftSize);
    const std::string rightContig = contigSeq.substr(endPos, rightSize);

    const int   searchWindow(500);
    const float mismatchRate(0.05f);
    const int   refAlignStart = align.beginPos;
    const int   refAlignEnd   = align.beginPos + apath_ref_length(apath);

    // search leftContig in the downstream of refStart
    const int         leftSearchStart     = std::max(0, refAlignEnd - searchWindow);
    const std::string refSeqForLeftSearch = refSeq.substr(leftSearchStart, (refAlignEnd - leftSearchStart));
    unsigned          occurrences = getQuerySeqMatchCount(refSeqForLeftSearch, leftContig, mismatchRate);

#ifdef DEBUG_CONTIG
    log_os << __FUNCTION__ << ": refSeqForLeftSearch: \n" << refSeqForLeftSearch << "\n";
    log_os << __FUNCTION__ << ": left contig has size " << leftSize << ":\n" << leftContig << "\n";
    log_os << __FUNCTION__ << ": left contig occurrences " << occurrences << "\n";
#endif
    if (occurrences > 1) return false;

    // search rightContig in the upstream of refEnd
    const int         rightSearchSize      = std::min(searchWindow, int(refSeq.length() - refAlignStart));
    const std::string refSeqForRightSearch = refSeq.substr(refAlignStart, rightSearchSize);
    occurrences = getQuerySeqMatchCount(refSeqForRightSearch, rightContig, mismatchRate);

#ifdef DEBUG_CONTIG
    log_os << __FUNCTION__ << ": refSeqForRightSearch: \n" << refSeqForRightSearch << "\n";
    log_os << __FUNCTION__ << ": right contig has size " << rightSize << ":\n" << rightContig << "\n";
    log_os << __FUNCTION__ << ": right contig occurrences " << occurrences << "\n";
#endif
    if (occurrences > 1) return false;
  }

  // check if an indel of at least minimum size is identified
  typedef std::pair<unsigned, unsigned> segment_t;
  std::vector<segment_t>                tmpseg(candidateSegments);
  candidateSegments.clear();
  for (const segment_t& segment : tmpseg) {
    for (unsigned i(segment.first); i <= segment.second; ++i) {
      if (((apath[i].type == INSERT) && (apath[i].length >= minCandidateIndelSize)) ||
          ((apath[i].type == DELETE) && (apath[i].length >= minCandidateIndelSize))) {
        candidateSegments.push_back(segment);
        break;
      }
    }
  }

  return (!candidateSegments.empty());
}

/// If a large insertion is not complete assembled, it must be assembled at least this far into either side
static const unsigned minSemiLargeInsertionLength(40);

/// \param[in] trimInsertLength remove extra length from the end of the contig
/// for the purpose of determining if the "unaligned" end is long enough
///
/// \return true if this is a left->right insert candidate
///
static bool isLargeInsertSegment(
    const AlignerBase<int>&  aligner,
    const ALIGNPATH::path_t& apath,
    unsigned&                contigOffset,
    unsigned&                refOffset,
    int&                     score,
    const unsigned           trimInsertLength = 0)
{
  using namespace ALIGNPATH;

  // Min length of aligned portion of contig
  static const unsigned minAlignReadLength(40);

  // Min length of unaligned portion of contig
  static const unsigned minExtendedReadLength(minSemiLargeInsertionLength);

  // Min reference length for alignment
  static const unsigned minAlignRefSpan(40);

  // Min fraction of optimal score in each contig sub-alignment
  static const float minScoreFrac(0.75);

  const unsigned pathSize(apath_read_length(apath));

  // first evaluate in the forward direction
  score = (std::max(0, getMaxPathScore(aligner.getScores(), apath, contigOffset, refOffset)));

#ifdef DEBUG_REFINER
  log_os << __FUNCTION__ << ": apath " << apath << "\n";
  log_os << __FUNCTION__ << ": score/ref/contig " << score << " " << refOffset << " " << contigOffset << "\n";
#endif

  if (refOffset < minAlignRefSpan) return false;
  if (contigOffset < minAlignReadLength) return false;

  assert(contigOffset <= pathSize);
  if ((pathSize - contigOffset) < (minExtendedReadLength + trimInsertLength)) return false;

  const int optimalScore(contigOffset * aligner.getScores().match);

  const float scoreFrac(static_cast<float>(score) / static_cast<float>(optimalScore));
  if (scoreFrac < minScoreFrac) return false;

  return true;
}

/// \return true if there is a large insert candidate
///
static bool isLargeInsertAlignment(
    const GlobalAligner<int> aligner, const ALIGNPATH::path_t& apath, LargeInsertionInfo& insertInfo)
{
  using namespace ALIGNPATH;

  insertInfo.isLeftCandidate =
      isLargeInsertSegment(aligner, apath, insertInfo.contigOffset, insertInfo.refOffset, insertInfo.score);

  if (insertInfo.isLeftCandidate) {
    return true;
  }

  ALIGNPATH::path_t apath_rev(apath);
  std::reverse(apath_rev.begin(), apath_rev.end());

  insertInfo.isRightCandidate = isLargeInsertSegment(
      aligner, apath_rev, insertInfo.contigOffset, insertInfo.refOffset, insertInfo.score);

  if (insertInfo.isRightCandidate) {
    const unsigned contigSize(apath_read_length(apath));
    const unsigned refSize(apath_ref_length(apath));
    insertInfo.contigOffset = contigSize - insertInfo.contigOffset;
    insertInfo.refOffset    = refSize - insertInfo.refOffset;
    return true;
  }

  return false;
}

/// \return true if the alignment represents an acceptable complete insertion:
///
static bool isFinishedLargeInsertAlignment(
    const GlobalAligner<int>             aligner,
    const ALIGNPATH::path_t&             apath,
    const std::pair<unsigned, unsigned>& insertSegment,
    const unsigned                       middleSize)
{
  using namespace ALIGNPATH;

  const path_t apath_left(apath.begin(), apath.begin() + insertSegment.second + 1);
  ;

  LargeInsertionInfo insertInfo;
  insertInfo.isLeftCandidate = isLargeInsertSegment(
      aligner, apath_left, insertInfo.contigOffset, insertInfo.refOffset, insertInfo.score, middleSize);

  path_t apath_rev(apath.begin() + insertSegment.first, apath.end());
  ;
  std::reverse(apath_rev.begin(), apath_rev.end());

  insertInfo.isRightCandidate = isLargeInsertSegment(
      aligner, apath_rev, insertInfo.contigOffset, insertInfo.refOffset, insertInfo.score, middleSize);

  return (insertInfo.isLeftCandidate && insertInfo.isRightCandidate);
}

/// Get the range over which an alignment element can vary with equal edit distance
///
/// \param[in] refRange range of the event (ie indel) of interest in reference coordinates
/// \param[in] readRange range of the event (ie indel) of interest in read coordinates
///
/// range coordinates are zero indexed and start at the first affected positions (so are not like vcf
/// coordinates) for instance:
////  the deletion 10M1D10M would have refRange(10,11), readRange(10,10)
////  the insertion 10M1I10M would have refRange(10,10), readRange(10,11)
///
static known_pos_range2 getVariantRange(
    const std::string&      ref,
    const known_pos_range2& refRange,
    const std::string&      read,
    const known_pos_range2& readRange)
{
#ifdef DEBUG_VARR
  log_os << __FUNCTION__ << ": refRange " << refRange << "\n";
  log_os << __FUNCTION__ << ": ref:\n";
  printSeq(ref, log_os);
  log_os << "\n";
  log_os << __FUNCTION__ << ": readRange " << readRange << "\n";
  log_os << __FUNCTION__ << ": read:\n";
  printSeq(read, log_os);
  log_os << "\n";
#endif

  // check how far we can slide to the right:
  const pos_t maxRightOffset(std::min(ref.size() - refRange.end_pos(), read.size() - readRange.end_pos()));
  pos_t       rightOffset(0);
  for (; rightOffset < maxRightOffset; ++rightOffset) {
    const char refSym(ref[refRange.begin_pos() + rightOffset]);
    const char readSym(read[readRange.begin_pos() + rightOffset]);
    if (refSym != readSym) break;
  }

  // check how far we can slide to the left:
  const pos_t minLeftOffset(std::max(-refRange.begin_pos(), -readRange.begin_pos()));
  pos_t       leftOffset(0);
  for (; leftOffset >= minLeftOffset; --leftOffset) {
    const char refSym(ref[refRange.end_pos() + leftOffset - 1]);
    const char readSym(read[readRange.end_pos() + leftOffset - 1]);
    if (refSym != readSym) break;
  }

#ifdef DEBUG_VARR
  log_os << __FUNCTION__ << ": left/right offset " << leftOffset << "/" << rightOffset << "\n";
#endif

  return known_pos_range2(leftOffset, rightOffset);
}

/// process smallSV alignment section into a usable sv candidate
static void setSmallCandSV(
    const reference_contig_segment&      ref,
    const std::string&                   contig,
    const Alignment&                     align,
    const std::pair<unsigned, unsigned>& segRange,
    SVCandidate&                         sv,
    const GSCOptions&                    opt)
{
#ifdef DEBUG_VARR
  log_os << __FUNCTION__ << ": align " << align << "\n";
  log_os << __FUNCTION__ << ": segRange [" << segRange.first << "," << segRange.second << "]\n";
  log_os << __FUNCTION__ << ": inputSV " << sv << "\n";
#endif
  sv.setPrecise();

  // get readRange and refRange, which are translations of segRange into
  // read and reference offsets:
  known_pos_range2 readRange;
  known_pos_range2 refRange;
  {
    using namespace ALIGNPATH;

    pos_t readPos(0);
    pos_t refPos(align.beginPos);

    const path_t&  apath(align.apath);
    const unsigned as(apath.size());
    for (unsigned i(0); i < as; ++i) {
      const path_segment& ps(apath[i]);
      if (i == segRange.first) {
        refRange.set_begin_pos(refPos);
        readRange.set_begin_pos(readPos);
      }

      if (is_segment_type_ref_length(ps.type)) refPos += ps.length;
      if (is_segment_type_read_length(ps.type)) readPos += ps.length;

      if (i == segRange.second) {
        refRange.set_end_pos(refPos);
        readRange.set_end_pos(readPos);
      }
    }
  }

  // by how many positions can the alignment position vary with the same alignment score?:
  const known_pos_range2 cipos(getVariantRange(ref.seq(), refRange, contig, readRange));

  // cipos for a precise variant is expected to start from 0 and extend forward zero to many bases
  if (cipos.begin_pos() != 0) {
    using namespace illumina::common;

    std::ostringstream oss;
    oss << "Attempting to convert alignment to sv candidate."
        << " contigSize: " << contig.size() << " alignment: " << align << " segments: [" << segRange.first
        << "," << segRange.second << "]\n"
        << "\treadRange: " << readRange << "\n"
        << "\trefRange: " << refRange << "\n"
        << "\tcipos: " << cipos << "\n";
    BOOST_THROW_EXCEPTION(GeneralException(oss.str()));
  }

  sv.bp1.state = SVBreakendState::RIGHT_OPEN;
  const pos_t beginPos(ref.get_offset() + refRange.begin_pos() - 1);
  sv.bp1.interval.range.set_range(beginPos, beginPos + cipos.end_pos() + 1);

  sv.bp2.state = SVBreakendState::LEFT_OPEN;
  const pos_t endPos(ref.get_offset() + refRange.end_pos());
  sv.bp2.interval.range.set_range(endPos, endPos + cipos.end_pos() + 1);
  sv.bp2.interval.tid = sv.bp1.interval.tid;

  sv.insertSeq = contig.substr(readRange.begin_pos(), readRange.size());

  // add CIGAR for all indels:
  sv.insertAlignment =
      ALIGNPATH::path_t(align.apath.begin() + segRange.first, align.apath.begin() + segRange.second + 1);

  // set assembled contig sequence, if option is specified
  if (opt.isOutputContig) {
    sv.contigSeq = contig;
  }
}

static known_pos_range2 getInsertTrim(
    const ALIGNPATH::path_t& apath, const std::pair<unsigned, unsigned>& segRange)
{
  assert(segRange.first <= segRange.second);

  using namespace ALIGNPATH;

  known_pos_range2 range;

  pos_t readPos(0);

  const unsigned as(apath.size());
  for (unsigned i(0); i < as; ++i) {
    const path_segment& ps(apath[i]);
    if (i == segRange.first) {
      range.set_begin_pos(readPos);
    }

    if (is_segment_type_read_length(ps.type)) readPos += ps.length;

    if (i == segRange.second) {
      range.set_end_pos(readPos);
      return range;
    }
  }

  assert(false && "segRange not found");
  return range;
}

/// Search for combinations of left and right-side insertion candidates to find a good insertion pair
static void processLargeInsertion(
    const SVCandidate&           sv,
    const pos_t                  leadingCut,
    const pos_t                  trailingCut,
    const GlobalAligner<int>&    largeInsertCompleteAligner,
    const std::vector<unsigned>& largeInsertionCandidateIndex,
    const std::set<pos_t>&       excludedPos,
    SVCandidateAssemblyData&     assemblyData,
    const GSCOptions&            opt)
{
  if (largeInsertionCandidateIndex.empty()) return;

#ifdef DEBUG_REFINER
  static const std::string logtag("processLargeInsertion: ");
  log_os << logtag << "starting large insertion search\n";
#endif

  bool     isLargeInsertionPair(false);
  unsigned largeInsertionLeftIndex(0);
  unsigned largeInsertionRightIndex(0);
  int      bestBreakDist(0);
  int      bestBreakScore(0);

  // try to pair up a large insertion candidate
  //
  // just do a dumb, all against all evaluation for now, if there's more than one left-right candidate set,
  // resolve according to (1) min ref distance and (2) best combined score
  static const int maxBreakDist(35);

  const unsigned candCount(largeInsertionCandidateIndex.size());
  for (unsigned candCount1(0); (candCount1 + 1) < candCount; ++candCount1) {
    const unsigned            candIndex1(largeInsertionCandidateIndex[candCount1]);
    const Alignment&          align1(assemblyData.smallSVAlignments[candIndex1].align);
    const LargeInsertionInfo& insert1(assemblyData.largeInsertInfo[candIndex1]);
    for (unsigned candCount2(candCount1 + 1); candCount2 < candCount; ++candCount2) {
      const unsigned            candIndex2(largeInsertionCandidateIndex[candCount2]);
      const Alignment&          align2(assemblyData.smallSVAlignments[candIndex2].align);
      const LargeInsertionInfo& insert2(assemblyData.largeInsertInfo[candIndex2]);
      if (!((insert1.isLeftCandidate && insert2.isRightCandidate) ||
            (insert2.isLeftCandidate && insert1.isRightCandidate)))
        continue;

      const int breakDist(std::abs(
          (long int)(align1.beginPos + insert1.refOffset) - (long int)(align2.beginPos + insert2.refOffset)));

      if (breakDist > maxBreakDist) continue;

      const int breakScore(insert1.score + insert2.score);

      const bool isBetterCandidate(
          (breakDist < bestBreakDist) || ((breakDist == bestBreakDist) && ((breakScore > bestBreakScore))));
      if ((!isLargeInsertionPair) || isBetterCandidate) {
        /// set new large insertion candidate:
        isLargeInsertionPair     = true;
        largeInsertionLeftIndex  = candIndex1;
        largeInsertionRightIndex = candIndex2;
        if (insert1.isRightCandidate) {
          std::swap(largeInsertionLeftIndex, largeInsertionRightIndex);
        }
        bestBreakDist  = breakDist;
        bestBreakScore = breakScore;

#ifdef DEBUG_REFINER
        log_os << logtag << "setting new large insertion candidate. Score: " << breakScore << "\n";
#endif
      }
    }
  }

  // no large insertion found:
  if (!isLargeInsertionPair) return;

  // found large insertion, insert this into data structures for downstream scoring/reporting:
  {
    const std::string& align1RefStr(assemblyData.bp1ref.seq());
    const unsigned     contigCount(assemblyData.contigs.size());

    static const std::string middle(
        "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN");
    const unsigned middleSize(middle.size());

    assemblyData.contigs.resize(contigCount + 1);
    assemblyData.smallSVAlignments.resize(contigCount + 1);
    assemblyData.smallSVSegments.resize(contigCount + 1);
    assemblyData.extendedContigs.resize(contigCount + 1);

    AssembledContig&                                   fakeContig(assemblyData.contigs[contigCount]);
    SVCandidateAssemblyData::SmallAlignmentResultType& fakeAlignment(
        assemblyData.smallSVAlignments[contigCount]);
    std::vector<std::pair<unsigned, unsigned>>& fakeSegments(assemblyData.smallSVSegments[contigCount]);
    std::string&                                fakeExtendedContig(assemblyData.extendedContigs[contigCount]);

    const AssembledContig& leftContig(assemblyData.contigs[largeInsertionLeftIndex]);
    const AssembledContig& rightContig(assemblyData.contigs[largeInsertionRightIndex]);

    fakeContig = leftContig;
    fakeContig.seq += (middle + rightContig.seq);

    const AssembledContig& constFakeContig(fakeContig);

    largeInsertCompleteAligner.align(
        constFakeContig.seq.begin(),
        constFakeContig.seq.end(),
        align1RefStr.begin() + leadingCut,
        align1RefStr.end() - trailingCut,
        fakeAlignment);

    fakeAlignment.align.beginPos += leadingCut;

#ifdef DEBUG_REFINER
    log_os << logtag << "large insertion fake alignment: " << fakeAlignment << "\n";
#endif

    fakeSegments.clear();
    getLargestInsertSegment(fakeAlignment.align.apath, middleSize, fakeSegments);

    // QC segments
    if ((1 != fakeSegments.size()) || (fakeSegments[0].second < fakeSegments[0].first)) {
      return;
    }

    // QC the resulting alignment:
    if (!isFinishedLargeInsertAlignment(
            largeInsertCompleteAligner, fakeAlignment.align.apath, fakeSegments[0], middleSize)) {
      return;
    }

    // final prep step: check left and right partial insert sequences -- this is a last chance to QC for
    // anomalies and get out:
    //
    const known_pos_range2 insertTrim(getInsertTrim(fakeAlignment.align.apath, fakeSegments[0]));
    {
      static const int minFlankSize(minSemiLargeInsertionLength);
      if ((insertTrim.begin_pos() + minFlankSize) > static_cast<pos_t>(leftContig.seq.size())) {
        return;
      }

      const pos_t rightOffset(leftContig.seq.size() + middle.size());
      if ((rightOffset + minFlankSize) > insertTrim.end_pos()) {
        return;
      }
    }

#ifdef DEBUG_REFINER
    log_os << logtag << "large insertion passed QC\n";
#endif

    getExtendedContig(fakeAlignment, fakeContig.seq, align1RefStr, fakeExtendedContig);

    // this section mostly imitates the regular SV build below, now that we've constructed our fake
    // contig/alignment
    SVCandidate newSV(sv);
    newSV.assemblyAlignIndex   = contigCount;
    newSV.assemblySegmentIndex = 0;
    setSmallCandSV(assemblyData.bp1ref, fakeContig.seq, fakeAlignment.align, fakeSegments[0], newSV, opt);

    /// check if this matches a fully assembled insertion:
    const pos_t startPos(newSV.bp1.interval.range.begin_pos());
    if (excludedPos.count(startPos)) return;

    newSV.isUnknownSizeInsertion = true;

    // final step: get left and right partial insert sequences:
    //
    assert(insertTrim.begin_pos() < static_cast<pos_t>(leftContig.seq.size()));
    newSV.unknownSizeInsertionLeftSeq = leftContig.seq.substr(insertTrim.begin_pos());

    const pos_t rightOffset(leftContig.seq.size() + middle.size());
    assert(rightOffset < insertTrim.end_pos());

    newSV.unknownSizeInsertionRightSeq = rightContig.seq.substr(0, (insertTrim.end_pos() - rightOffset));

    assemblyData.svs.push_back(newSV);
  }
}

SVCandidateAssemblyRefiner::SVCandidateAssemblyRefiner(
    const GSCOptions&                   opt,
    const bam_header_info&              header,
    const AllSampleReadCounts&          counts,
    std::shared_ptr<EdgeRuntimeTracker> edgeTrackerPtr)
  : _opt(opt),
    _header(header),
    _smallSVAssembler(
        opt.scanOpt,
        opt.refineOpt.smallSVAssembleOpt,
        opt.alignFileOpt,
        opt.referenceFilename,
        opt.statsFilename,
        opt.chromDepthFilename,
        header,
        counts,
        opt.isRNA,
        edgeTrackerPtr->remoteReadRetrievalTime),
    _spanningAssembler(
        opt.scanOpt,
        (opt.isRNA ? opt.refineOpt.RNAspanningAssembleOpt : opt.refineOpt.spanningAssembleOpt),
        opt.alignFileOpt,
        opt.referenceFilename,
        opt.statsFilename,
        opt.chromDepthFilename,
        header,
        counts,
        opt.isRNA,
        edgeTrackerPtr->remoteReadRetrievalTime),
    _largeSVAligner(opt.refineOpt.largeSVAlignScores, opt.refineOpt.largeGapOpenScore),
    _largeInsertEdgeAligner(opt.refineOpt.largeInsertEdgeAlignScores),
    _largeInsertCompleteAligner(opt.refineOpt.largeInsertCompleteAlignScores),
    _spanningAligner(opt.refineOpt.spanningAlignScores, opt.refineOpt.jumpScore),
    _RNASpanningAligner(
        opt.refineOpt.RNAspanningAlignScores,
        opt.refineOpt.RNAJumpScore,
        opt.refineOpt.RNAIntronOpenScore,
        opt.refineOpt.RNAIntronOffEdgeScore),
    _contigFilterAlignmentScores(opt.refineOpt.contigFilterScores)
{
}

void SVCandidateAssemblyRefiner::getCandidateAssemblyData(
    const SVCandidate& sv, const bool isFindLargeInsertions, SVCandidateAssemblyData& assemblyData) const
{
#ifdef DEBUG_REFINER
  log_os << __FUNCTION__ << ": START sv " << sv;
#endif

  assemblyData.clear();

  // separate the problem into different assembly categories:
  //
  if (isSpanningSV(sv)) {
    // record the spanning status of the original low-resolution candidate:
    assemblyData.isCandidateSpanning = true;

    // this case assumes two suspected breakends with a direction to each, most common large scale SV case:
    getJumpAssembly(sv, isFindLargeInsertions, assemblyData);
  } else if (isComplexSV(sv)) {
    // record the spanning status of the original low-resolution candidate:
    assemblyData.isCandidateSpanning = false;

    // this case assumes a single-interval local assembly, this is the most common case for small-scale
    // SVs/indels
    getSmallSVAssembly(sv, isFindLargeInsertions, assemblyData);
  } else {
    log_os << "Unknown candidate SV: " << sv << "\n";
    assert(false && "Unknown candidate SV type");
  }
}

/// Represents a stretch of reference sequence excluded from alignment by the kmer matcher.
struct exclusion_block {
  exclusion_block(const unsigned s, const unsigned l, const unsigned sp) : start(s), length(l), nSpacer(sp) {}
  unsigned start;    // Start of excluded region
  unsigned length;   // Nunmber of bp excluded
  unsigned nSpacer;  // Number of 'N' added in place of excluded sequence
};

/// Translate a reduced reference position to the original reference coordinates
static unsigned translateMaskedPos(const std::vector<exclusion_block>& exclBlocks, const unsigned maskedPos)
{
  int offset = 0;
  for (const auto& cblock : exclBlocks) {
    if (cblock.start > (offset + maskedPos)) break;
    offset += cblock.length - cblock.nSpacer;
  }
  return offset + maskedPos;
}

/// Translate an alignment made against a reduced reference to the original reference coordinates
static bool translateMaskedAlignment(Alignment& align, const std::vector<exclusion_block>& exclBlocks)
{
  using namespace ALIGNPATH;
#ifdef DEBUG_KMER
  log_os << __FUNCTION__ << " original: " << align << "\n";
#endif
  path_t newPath;
  pos_t  cpos = align.beginPos;
  for (const path_segment& seg : align.apath) {
    if (!is_segment_type_ref_length(seg.type)) {
      newPath.push_back(seg);
    } else {
      const unsigned length =
          translateMaskedPos(exclBlocks, cpos + seg.length) - translateMaskedPos(exclBlocks, cpos);
      if (is_segment_align_match(seg.type) && (length != seg.length)) return false;
#ifdef DEBUG_KMER
      log_os << __FUNCTION__ << " SEGMENT " << seg.type << " " << seg.length << "\n";
      log_os << "\tlength " << length << "\n";
      log_os << "\tcpos " << cpos << "\n";
#endif
      cpos += seg.length;
      newPath.emplace_back(seg.type, length);
    }
  }
  if (align.apath.size() > 0) {
    align.beginPos = translateMaskedPos(exclBlocks, align.beginPos);
    align.apath    = newPath;
  }
#ifdef DEBUG_KMER
  log_os << __FUNCTION__ << " final: " << align << "\n";
#endif
  return true;
}

namespace {

/// Returns a reduced reference sequence where long stretches without kmer matches to the contig are removed
template <typename SymIter>
std::string kmerMaskReference(
    const SymIter                 refSeqStart,
    const SymIter                 refSeqEnd,
    const std::string&            contig,
    const int                     nSpacer,
    std::vector<exclusion_block>& exclBlocks)
{
  // Hash all kmers in the contig
  static const int                merSize(10);
  std::unordered_set<std::string> contigHash;
  for (unsigned contigMerIndex(0); contigMerIndex < (contig.size() - (merSize - 1)); ++contigMerIndex) {
    contigHash.insert(contig.substr(contigMerIndex, merSize));
  }
  // Mask the reference (and keep track of excluded regions for coordinate translation later)
  static const int minExclusion(1000);
  static const int padding(50);  // Amount of sequence included around each kmer hit.
  std::string      maskedRef;
  const SymIter    maxRef       = refSeqEnd - (merSize - 1);
  SymIter          potExclStart = refSeqStart;
  SymIter          inclStart    = refSeqStart;
  for (SymIter refIt = refSeqStart; refIt != maxRef; refIt++) {
    if (contigHash.count(std::string(refIt, refIt + merSize)) != 0) {
      if ((refIt - potExclStart) > (minExclusion + padding)) {
        unsigned spacer(0);
        if (potExclStart > refSeqStart) {
          maskedRef.append(std::string(inclStart, potExclStart));
          maskedRef.append(nSpacer, 'N');
          spacer = nSpacer;
        }
        inclStart = refIt - padding;
        exclBlocks.emplace_back(potExclStart - refSeqStart, inclStart - potExclStart, spacer);
      }
      potExclStart = refIt + padding;
    }
  }
  maskedRef.append(std::string(inclStart, std::min(maxRef, potExclStart)));
#ifdef DEBUG_KMER
  log_os << __FUNCTION__ << "Reduced to " << maskedRef << '\n';
  log_os << __FUNCTION__ << " exclBlocks\n\t";
  for (const auto block : exclBlocks)
    log_os << " " << block.start << ":" << block.length << ":" << block.nSpacer;
  log_os << "\n";
#endif
  if (maskedRef.empty()) maskedRef.append(nSpacer, 'N');
  return maskedRef;
}

}  // namespace

/// Convert jump alignment results into an SVCandidate
///
static void generateRefinedSVCandidateFromJumpAlignment(
    const SVCandidateAssemblyData& assemblyData, SVCandidate& sv)
{
  const SVCandidateAssemblyData::JumpAlignmentResultType& align(
      assemblyData.spanningAlignments[assemblyData.bestAlignmentIndex]);

  // first get each alignment associated with the correct breakend:
  const Alignment* bp1AlignPtr(&align.align1);
  const Alignment* bp2AlignPtr(&align.align2);

  if (assemblyData.bporient.isBp2AlignedFirst) std::swap(bp1AlignPtr, bp2AlignPtr);

  // summarize usable output information in a second SVBreakend
  // object -- this is the 'refined' sv:
  sv.assemblyAlignIndex   = assemblyData.bestAlignmentIndex;
  sv.assemblySegmentIndex = 0;

  sv.setPrecise();

  adjustAssembledBreakend(
      *bp1AlignPtr,
      (!assemblyData.bporient.isBp2AlignedFirst),
      align.jumpRange,
      assemblyData.bp1ref,
      assemblyData.bporient.isBp1Reversed,
      sv.bp1);
  adjustAssembledBreakend(
      *bp2AlignPtr,
      (assemblyData.bporient.isBp2AlignedFirst),
      align.jumpRange,
      assemblyData.bp2ref,
      assemblyData.bporient.isBp2Reversed,
      sv.bp2);
}

/// Convert jump alignment results into an SVCandidate and add all
/// extra data required for VCF output
///
static void generateRefinedVCFSVCandidateFromJumpAlignment(
    const SVCandidateAssemblyData& assemblyData, SVCandidate& sv, const GSCOptions& opt)
{
  generateRefinedSVCandidateFromJumpAlignment(assemblyData, sv);

  const AssembledContig& contig(assemblyData.contigs[assemblyData.bestAlignmentIndex]);
  const SVCandidateAssemblyData::JumpAlignmentResultType& align(
      assemblyData.spanningAlignments[assemblyData.bestAlignmentIndex]);

  // fill in insertSeq:
  sv.insertSeq.clear();
  if (align.jumpInsertSize > 0) {
    getFwdStrandInsertSegment(align, contig.seq, assemblyData.bporient.isBp1Reversed, sv.insertSeq);
  }

  // fill in contigSeq, only when "--outputConfig" is specified
  if (opt.isOutputContig) {
    sv.contigSeq = contig.seq;
  }

  // add CIGAR for any simple (insert/delete) cases:
  addCigarToSpanningAlignment(sv);
}

/// \brief Helper function for isJumpAlignmentQCFail
/// \return True if the jump alignment segment fails basic QC checks
static bool isJumpAlignmentSegmentQCFail(const Alignment& alignment)
{
  static const unsigned minAlignRefSpan(20);
  return ((!alignment.isAligned()) || (apath_ref_length(alignment.apath) < minAlignRefSpan));
}

/// \brief A fast QC check of the jump alignment to ensure it spans the two breakend locations
///
/// This jump alignment routine does a less thorough, but simpler assesement of the jump alignment intended
/// for a quick QC screen
///
/// \return True if the jump alignment fails basic QC checks
static bool isJumpAlignmentQCFail(const JumpAlignmentResult<int>& jumpAlignment)
{
  return (
      isJumpAlignmentSegmentQCFail(jumpAlignment.align1) ||
      isJumpAlignmentSegmentQCFail(jumpAlignment.align2));
}

/// \brief Check if the jump alignment result is low quality
///
/// Check the min size and fraction of optimal alignment score for each of the two sub-alignments on each side
/// of the breakpoint
///
/// Note this is done for multiple values -- the lower value is motivated by cases where a second breakpoint
/// exists near to the target breakpoint -- the higher value is motivated by cases with some alignment
/// 'messiness' near the breakpoint that stabilizes as we move farther away
/// TODO change iterative refspan to a single consistent alignment criteria
///
/// \param[in] alignmentScores Scores used for the purpose of assessing the quality of the jump aligner's two
/// subalignments on each side of the breakend
///
/// \return True if the jump alignment result is low quality and should be filtered out
static bool isLowQualityJumpAlignment(
    const JumpAlignmentResult<int>& alignment, const AlignmentScores<int>& alignmentScores, const bool isRNA)
{
  using namespace ALIGNPATH;
  bool                  isLowQualityAlign1(true);
  bool                  isLowQualityAlign2(true);
  static const unsigned spanSet[] = {75, 100, 200};
  // Use shorter regions for RNA contig alignment QC since reads and inserts
  // (and hence assembled contigs) are typically shorter
  const unsigned spanSetRna[] = {36, 75, 100};

  for (const unsigned maxQCRefSpan : isRNA ? spanSetRna : spanSet) {
    const unsigned qcSpan1 = maxQCRefSpan + (isRNA ? apath_spliced_length(alignment.align1.apath) : 0);
    if (!isLowQualitySpanningSVAlignment(qcSpan1, alignmentScores, true, isRNA, alignment.align1.apath)) {
      isLowQualityAlign1 = false;
    }
    const unsigned qcSpan2 = maxQCRefSpan + (isRNA ? apath_spliced_length(alignment.align2.apath) : 0);
    if (!isLowQualitySpanningSVAlignment(qcSpan2, alignmentScores, false, isRNA, alignment.align2.apath)) {
      isLowQualityAlign2 = false;
    }
  }
  return (isLowQualityAlign1 || isLowQualityAlign2);
}

/// Filter fusion contigs and select the 'best' one, based on alignment score and supporting read count
static bool selectJumpContigRNA(
    SVCandidateAssemblyData& assemblyData, const AlignmentScores<int>& alignmentScores)
{
  static const bool isRNA(true);
  const unsigned    contigCount(assemblyData.contigs.size());

  std::vector<int> goodContigIndices;
  for (unsigned contigIndex = 0; contigIndex < contigCount; contigIndex++) {
    const JumpAlignmentResult<int>& alignment(assemblyData.spanningAlignments[contigIndex]);
#ifdef DEBUG_REFINER
    log_os << __FUNCTION__ << ": Checking contig alignment: " << contigIndex << "\n";
#endif
    if (isJumpAlignmentQCFail(alignment)) continue;
#ifdef DEBUG_REFINER
    log_os << __FUNCTION__ << ": contig alignment initial okay: " << contigIndex << "\n";
#endif
    if (isLowQualityJumpAlignment(alignment, alignmentScores, isRNA)) continue;
#ifdef DEBUG_REFINER
    log_os << __FUNCTION__ << ": contig alignment okay: " << contigIndex << "\n";
#endif
    goodContigIndices.push_back(contigIndex);
  }
  if (goodContigIndices.empty()) return false;
  // Find the highest alignment score
  int      maxAlnScore = 0;
  unsigned selectedContigIndex(goodContigIndices.front());
  for (unsigned index : goodContigIndices) {
    if (assemblyData.spanningAlignments[index].score > maxAlnScore) {
      maxAlnScore         = assemblyData.spanningAlignments[index].score;
      selectedContigIndex = index;
    }
  }
  // Pick the contig with the most supporting reads that has an alignment score at least half as high as the
  // highest-scoring contig
  for (unsigned index : goodContigIndices) {
    const bool sufficientScore(assemblyData.spanningAlignments[index].score * 2 > maxAlnScore);
    const bool moreReads(
        assemblyData.contigs[index].supportReads.size() >
        assemblyData.contigs[selectedContigIndex].supportReads.size());
    if (sufficientScore && moreReads) selectedContigIndex = index;
  }
#ifdef DEBUG_REFINER
  log_os << __FUNCTION__ << ": selected contig: " << selectedContigIndex << "\n";
#endif
  assemblyData.bestAlignmentIndex = selectedContigIndex;
  return true;
}

/// Filter breakpoint contigs for large SV candidates
/// Select the 'best' one based on alignment score and check alignment on selected contig
// TODO Consider making this more like the RNA case, e.g. all alignment checks before selection and pick the
// best passing one.
static bool selectJumpContigDNA(
    SVCandidateAssemblyData& assemblyData, const AlignmentScores<int>& alignmentScores)
{
  static const bool isRNA(false);
  const auto&       alignments(assemblyData.spanningAlignments);
  const unsigned    contigCount(assemblyData.contigs.size());

  int maxAlignContigIndex(-1);
  for (unsigned contigIndex = 0; contigIndex < contigCount; contigIndex++) {
    const JumpAlignmentResult<int>& alignment(alignments[contigIndex]);
#ifdef DEBUG_REFINER
    log_os << __FUNCTION__ << ": Checking contig alignment: " << contigIndex << "\n";
#endif
    if (isJumpAlignmentQCFail(alignment)) continue;
#ifdef DEBUG_REFINER
    log_os << __FUNCTION__ << ": contig alignment initial okay: " << contigIndex << "\n";
#endif
    // Find the contig with the highest alignment score
    if ((maxAlignContigIndex == -1) ||
        (alignments[contigIndex].score > alignments[maxAlignContigIndex].score)) {
      maxAlignContigIndex = contigIndex;
    }
  }
#ifdef DEBUG_REFINER
  log_os << __FUNCTION__ << ": selected contig: " << maxAlignContigIndex << "\n";
#endif
  if ((maxAlignContigIndex == -1) ||
      (isLowQualityJumpAlignment(alignments[maxAlignContigIndex], alignmentScores, isRNA))) {
    return false;
  }
  // ok, passed QC -- mark the high-scoring alignment as usable for
  // hypothesis refinement:
  assemblyData.bestAlignmentIndex = maxAlignContigIndex;
  return true;
}

/// Store auxilary alignment info
struct AlignData {
  unsigned bp1LeadingTrim  = 0;
  unsigned bp1TrailingTrim = 0;
  unsigned bp2LeadingTrim  = 0;
  unsigned bp2TrailingTrim = 0;

  pos_t align1LeadingCut  = 0;
  pos_t align1TrailingCut = 0;
  pos_t align2LeadingCut  = 0;
  pos_t align2TrailingCut = 0;
};

/// Assemble candidate contigs for large SV candidates
///
/// \param[in] alignData Initialized auxilary alignment info for sequence trimming
///
/// \param[in] assemblyData Assembly data for a SV candidate that will be updated to include contigs assembled
/// and reference sequences
///
/// \Return true if the reference region is valid and the assembly procedure is completed
///
bool static assembleJumpContigs(
    const GSCOptions&           opt,
    const SVCandidate&          sv,
    const pos_t                 extraRefEdgeSize,
    const pos_t                 extraRefSplitSize,
    const bam_header_info&      header,
    const SVCandidateAssembler& spanningAssembler,
    AlignData&                  alignData,
    SVCandidateAssemblyData&    assemblyData)
{
  assemblyData.isSpanning = true;
  BPOrientation& bporient(assemblyData.bporient);

  bporient.isBp1First              = sv.isForward();
  bporient.isTranscriptStrandKnown = sv.isTranscriptStrandKnown();
  if (opt.isRNA) {
    bporient.isBp1First =
        !sv.isForward();  // RNA-seq reads generate candidates in the opposite direction of the RNA
  }
  //
  // based on sv candidate, we classify the expected relationship
  // between the contig and the sv breakends:
  //
  if (sv.bp1.state != sv.bp2.state) {
    // if there's one right-open breakend and one left-open breakend, no matter the bp1/bp2 chromosome and
    // relative bp1/bp2 order etc. we:
    // 1. don't need to do any read/reference reversals
    // 2. always treat the right-open breakend as the first alignment region in order:
    //
    if (sv.bp2.state == SVBreakendState::RIGHT_OPEN) {
      bporient.isBp2AlignedFirst = true;
    }
  } else {
    // If both breakends open in the same direction, then:
    // 1. the reads from one breakend need to be reversed
    // 2. the reference from that same breakend needs to be reversed
    // 3. Treat the un-reversed RIGHT_OPEN or reversed LEFT_OPEN as the first alignment region in order
    //      Note that in the scheme below, we chose which bp to reverse so that no-reordering is required
    //
    if (sv.bp1.state == SVBreakendState::RIGHT_OPEN) {
      bporient.isBp2Reversed = true;
    } else {
      bporient.isBp1Reversed = true;
    }
  }

  // there's always a small chance that our region could fall
  // completely off the edge of the reference b/c of circular
  // chromosomes, this can't be treated as a bug -- it's a
  // legitimate breakend hypothesis that we just aren't setup to
  // handle correctly, so we punt this case:
  if (!isRefRegionValid(header, sv.bp1.interval)) return false;
  if (!isRefRegionValid(header, sv.bp2.interval)) return false;

  // next we extract the reference sequence around both breakends
  //
  // the 'trim' values below refer to the difference between the
  // breakend reference region requested and the region returned
  // after accounting for chromosome edges. The trim values will
  // almost always be zero for large chromosomes.
  //
  const pos_t extraRefSize(extraRefEdgeSize + extraRefSplitSize);
  getSVReferenceSegments(
      opt.referenceFilename,
      header,
      extraRefSize,
      sv,
      assemblyData.bp1ref,
      assemblyData.bp2ref,
      alignData.bp1LeadingTrim,
      alignData.bp1TrailingTrim,
      alignData.bp2LeadingTrim,
      alignData.bp2TrailingTrim);

  // The *Cut values below represent sequence which will be removed from the edges of the reference region for
  // each breakend. In most cases this will equal extraRefSplitSize. Sometimes these values are forced to be
  // shorter because we didn't retrieve as much reference sequence as targeted.
  alignData.align1LeadingCut = std::max(0, extraRefSplitSize - static_cast<pos_t>(alignData.bp1LeadingTrim));
  alignData.align1TrailingCut =
      std::max(0, extraRefSplitSize - static_cast<pos_t>(alignData.bp1TrailingTrim));
  alignData.align2LeadingCut = std::max(0, extraRefSplitSize - static_cast<pos_t>(alignData.bp2LeadingTrim));
  alignData.align2TrailingCut =
      std::max(0, extraRefSplitSize - static_cast<pos_t>(alignData.bp2TrailingTrim));

  // assemble contig(s) spanning the breakend:
  spanningAssembler.assembleSpanningSVCandidate(
      sv.bp1,
      sv.bp2,
      bporient.isBp1Reversed,
      bporient.isBp2Reversed,
      assemblyData.bp1ref,
      assemblyData.bp2ref,
      assemblyData.contigs);

  return true;
}

/// Align contigs of large SV candidates to reference
///
/// \param alignData Auxilary alignment info for sequence trimming initialized in contig assembly that will be
/// updated during contig alignment
///
void static alignJumpContigs(
    const GSCOptions&                   opt,
    const SVCandidate&                  sv,
    const GlobalJumpAligner<int>&       spanningAligner,
    const GlobalJumpIntronAligner<int>& RNASpanningAligner,
    AlignData&                          alignData,
    SVCandidateAssemblyData&            assemblyData)
{
  BPOrientation& bporient(assemblyData.bporient);
  std::string    bp1refSeq = assemblyData.bp1ref.seq();
  std::string    bp2refSeq = assemblyData.bp2ref.seq();
  if (bporient.isBp1Reversed) {
    reverseCompStr(bp1refSeq);
    std::swap(alignData.align1LeadingCut, alignData.align1TrailingCut);
  }
  if (bporient.isBp2Reversed) {
    reverseCompStr(bp2refSeq);
    std::swap(alignData.align2LeadingCut, alignData.align2TrailingCut);
  }

  const std::string* align1RefStrPtr = &bp1refSeq;
  const std::string* align2RefStrPtr = &bp2refSeq;
  if (bporient.isBp2AlignedFirst) {
    std::swap(align1RefStrPtr, align2RefStrPtr);

    std::swap(alignData.align1LeadingCut, alignData.align2LeadingCut);
    std::swap(alignData.align1TrailingCut, alignData.align2TrailingCut);
  }

#ifdef DEBUG_REFINER
  log_os << __FUNCTION__ << ": align1RefSize/Seq: " << align1RefStrPtr->size() << '\n';
  printSeq(*align1RefStrPtr, log_os);
  log_os << '\n';
  log_os << __FUNCTION__ << ": align2Refsize/Seq: " << align2RefStrPtr->size() << '\n';
  printSeq(*align2RefStrPtr, log_os);
  log_os << '\n';
#endif

  const unsigned contigCount(assemblyData.contigs.size());
  // make sure an alignment object exists for every contig, even if it's empty
  assemblyData.spanningAlignments.resize(contigCount);

#ifdef DEBUG_REFINER
  log_os << __FUNCTION__ << ": contigCount: " << contigCount << "\n";
  for (unsigned contigIndex(0); contigIndex < contigCount; ++contigIndex) {
    const AssembledContig& contig(assemblyData.contigs[contigIndex]);
    log_os << __FUNCTION__ << ": contigIndex: " << contigIndex << " contig: " << contig;
  }
#endif

  for (unsigned contigIndex(0); contigIndex < contigCount; ++contigIndex) {
    const AssembledContig& contig(assemblyData.contigs[contigIndex]);

#ifdef DEBUG_REFINER
    log_os << __FUNCTION__ << ": start aligning contigIndex: " << contigIndex << "\n";
#endif

    JumpAlignmentResult<int>& alignment(assemblyData.spanningAlignments[contigIndex]);

    if (opt.isRNA) {
#ifdef DEBUG_REFINER
      log_os << __FUNCTION__ << " RNA alignment\n";
#endif
      static const int             nSpacer(25);
      std::vector<exclusion_block> exclBlocks1;
      const std::string            cutRef1 = kmerMaskReference(
          align1RefStrPtr->begin() + alignData.align1LeadingCut,
          align1RefStrPtr->end() - alignData.align1TrailingCut,
          contig.seq,
          nSpacer,
          exclBlocks1);
      std::vector<exclusion_block> exclBlocks2;
      const std::string            cutRef2 = kmerMaskReference(
          align2RefStrPtr->begin() + alignData.align2LeadingCut,
          align2RefStrPtr->end() - alignData.align2TrailingCut,
          contig.seq,
          nSpacer,
          exclBlocks2);
#ifdef DEBUG_REFINER
      log_os << __FUNCTION__ << " Kmer-masked references\n";
      log_os << "\t ref Lengths " << align1RefStrPtr->size() << " " << align2RefStrPtr->size() << "\n";
      log_os << "\t cutref Lengths " << cutRef1.size() << " " << cutRef2.size() << "\n";
#endif

      bool bp1RnaStrandFw;  // Is the RNA fusion transcript on the forward strand at bp1
      bool bp2RnaStrandFw;
      if (bporient.isBp1First) {
        bp1RnaStrandFw = (sv.bp1.state == SVBreakendState::RIGHT_OPEN);
        bp2RnaStrandFw = (sv.bp2.state == SVBreakendState::LEFT_OPEN);
      } else {
        bp1RnaStrandFw = (sv.bp1.state == SVBreakendState::LEFT_OPEN);
        bp2RnaStrandFw = (sv.bp2.state == SVBreakendState::RIGHT_OPEN);
      }
      // Should we look for the splice motif on the fw or rev strand in the bp1 ref seq
      bool bp1Fw = (bporient.isBp1Reversed != bp1RnaStrandFw);
      bool bp2Fw = (bporient.isBp2Reversed != bp2RnaStrandFw);
      // bp1 and bp2 sequences have been swapped above
      if (bporient.isBp2AlignedFirst) std::swap(bp1Fw, bp2Fw);

#ifdef DEBUG_REFINER
      log_os << __FUNCTION__ << " isTranscriptStrandKnown: " << bporient.isTranscriptStrandKnown
             << "; bp1Fw: " << bp1Fw << " ; bp2Fw: " << bp2Fw << '\n';
#endif
      RNASpanningAligner.align(
          contig.seq.begin(),
          contig.seq.end(),
          cutRef1.begin(),
          cutRef1.end(),
          cutRef2.begin(),
          cutRef2.end(),
          bp1Fw,
          bp2Fw,
          bporient.isTranscriptStrandKnown,
          alignment);

#ifdef DEBUG_REFINER
      log_os << __FUNCTION__ << " Masked 1: " << alignment.align1 << '\n';
      log_os << __FUNCTION__ << " Masked 2: " << alignment.align2 << '\n';
#endif
      if (!(translateMaskedAlignment(alignment.align1, exclBlocks1) &&
            translateMaskedAlignment(alignment.align2, exclBlocks2))) {
#ifdef DEBUG_REFINER
        log_os << __FUNCTION__ << " Failed to fix kmer-masked alignment\n";
#endif
        alignment.align1.clear();
        alignment.align2.clear();
      }
#ifdef DEBUG_REFINER
      log_os << __FUNCTION__ << " Fixed 1: " << alignment.align1 << '\n';
      log_os << __FUNCTION__ << " Fixed 2: " << alignment.align2 << '\n';
#endif
    } else {
#ifdef DEBUG_REFINER
      log_os << __FUNCTION__ << " Ref1 for alignment: "
             << bp1refSeq.substr(
                    alignData.align1LeadingCut,
                    bp1refSeq.size() - alignData.align1LeadingCut - alignData.align1TrailingCut)
             << '\n';
      log_os << __FUNCTION__ << " Ref2 for alignment: "
             << bp2refSeq.substr(
                    alignData.align2LeadingCut,
                    bp2refSeq.size() - alignData.align2LeadingCut - alignData.align2TrailingCut)
             << '\n';
#endif
      spanningAligner.align(
          contig.seq.begin(),
          contig.seq.end(),
          align1RefStrPtr->begin() + alignData.align1LeadingCut,
          align1RefStrPtr->end() - alignData.align1TrailingCut,
          align2RefStrPtr->begin() + alignData.align2LeadingCut,
          align2RefStrPtr->end() - alignData.align2TrailingCut,
          alignment);

      const bool  hasJumpInsert(alignment.jumpInsertSize > 0);
      const pos_t minAlignBuffer(5);
      const pos_t ref1EndPos(
          align1RefStrPtr->size() - alignData.align1LeadingCut - alignData.align1TrailingCut - 1);
      const pos_t align1EndPos(alignment.align1.beginPos + apath_ref_length(alignment.align1.apath));
      // Breakend from align1 is close to the end position of reference sequence
      const bool isBp1CloseToRef1End(ref1EndPos - align1EndPos < minAlignBuffer);
      // Breakend from align2 is close to the start position of reference sequence
      const bool isBp2CloseToRef2Start(alignment.align2.beginPos < minAlignBuffer);
      // When there is any breakend close to either end of reference sequence,
      // and when the jump alignment includes an insert,
      // the contig need be realigned with longer reference sequences
      // (i.e. the reference sequences without cutting).
      //
      // Note that the reference sequences without cutting is expanded at most by extraRefSplitSize(default
      // 100bp). Therefore this one-time expansion may not resolve the very rare case, where the breakend is
      // located more than extraRefSplitSize away from the reference end.
      if (hasJumpInsert && (isBp1CloseToRef1End || isBp2CloseToRef2Start)) {
        alignData.align1LeadingCut  = 0;
        alignData.align1TrailingCut = 0;
        alignData.align2LeadingCut  = 0;
        alignData.align2TrailingCut = 0;
#ifdef DEBUG_REFINER
        log_os << __FUNCTION__ << " Aglignment: " << alignment << "\n"
               << "The breakend is close to the start/end position of reference sequence!\n";
        log_os << __FUNCTION__ << " Realign with longer reference sequences: \n"
               << " Ref1 for alignment: " << bp1refSeq << '\n'
               << " Ref2 for alignment: " << bp2refSeq << '\n';
#endif
        spanningAligner.align(
            contig.seq.begin(),
            contig.seq.end(),
            align1RefStrPtr->begin(),
            align1RefStrPtr->end(),
            align2RefStrPtr->begin(),
            align2RefStrPtr->end(),
            alignment);
      }
    }

    alignment.align1.beginPos += alignData.align1LeadingCut;
    alignment.align2.beginPos += alignData.align2LeadingCut;

    std::string extendedContig;
    getExtendedContig(alignment, contig.seq, *align1RefStrPtr, *align2RefStrPtr, extendedContig);
    assemblyData.extendedContigs.push_back(extendedContig);

#ifdef DEBUG_REFINER
    log_os << __FUNCTION__ << ": contigIndex: " << contigIndex << " alignment: " << alignment;
    if (alignment.align1.isAligned()) {
      std::string bp1Seq, bp2Seq, insertSeq;
      getFwdStrandQuerySegments(
          alignment,
          contig.seq,
          bporient.isBp2AlignedFirst,
          bporient.isBp1Reversed,
          bporient.isBp2Reversed,
          bp1Seq,
          bp2Seq,
          insertSeq);
      log_os << __FUNCTION__ << "\tbp1seq_fwd: " << bp1Seq << "\n";
      log_os << __FUNCTION__ << "\tinsseq_fwd: " << insertSeq << "\n";
      log_os << __FUNCTION__ << "\tbp2seq_fwd: " << bp2Seq << "\n";
    }
#endif
  }
}

void SVCandidateAssemblyRefiner::getJumpAssembly(
    const SVCandidate& sv, const bool isFindLargeInsertions, SVCandidateAssemblyData& assemblyData) const
{
#ifdef DEBUG_REFINER
  log_os << __FUNCTION__ << ": START\n";
  if (_opt.isRNA) {
    log_os << __FUNCTION__ << ": RNA\n";
  }
#endif

  // This determines by how much we extend the reference sequence
  // around the breakend region for all discovery and scoring
  // operations. It is possible to discover breakends and small
  // indels in this expanded region.
  //
  const pos_t extraRefEdgeSize(_opt.isRNA ? 25000 : 250);
  // This determines by how much we extend the reference sequence
  // around the breakend region for all operations except for the
  // initial alignment of the contig back to the reference.
  //
  // The primary motivation for this value is to improve our ability
  // to find reads which support the breakend when quality scoring
  // takes place (a subsequent step outside of this function), but
  // without expanding the regions where a breakend can possibly be
  // found (ie. without the risk of additional false positives)
  //
  // The extra reference sequence is used for:
  // - Extraction of breakend associated reads (reference sequence
  // used to find reads which poorly match the reference in this
  // case)
  // - Contig extension after alignment, this means that the
  // extended region will be used for read support scoring (later
  // during Q-value generation)
  //
  // The extra reference sequence
  // - is removed for the initial contig alignment, this means that
  // actual breakpoint discovery will not occur in the extended region
  // - will be added back for contig realignment if the initial
  // alignment leads to breakends close to reference ends.
  //
  const pos_t extraRefSplitSize(100);

  const pos_t extraRefSize(extraRefEdgeSize + extraRefSplitSize);

  // if the breakends have a simple insert/delete orientation and
  // the alignment regions overlap, then handle this case as a local
  // assembly problem:
  if (sv.bp1.interval.tid == sv.bp2.interval.tid) {
    if (!SVBreakendState::isSameOrientation(sv.bp1.state, sv.bp2.state)) {
      const SV_TYPE::index_t svType(getSVType(sv));
      if ((svType == SV_TYPE::INDEL) || (svType == SV_TYPE::COMPLEX)) {
        if (isRefRegionOverlap(_header, extraRefSize, sv)) {
#ifdef DEBUG_REFINER
          log_os << __FUNCTION__
                 << ": Candidate breakends regions are too close, transferring problem to local assembler\n";
#endif
          // transform SV into a single region format:
          SVCandidate singleSV = sv;
          singleSV.bp1.state   = SVBreakendState::COMPLEX;
          singleSV.bp2.state   = SVBreakendState::UNKNOWN;
          singleSV.bp1.interval.range.merge_range(sv.bp2.interval.range);

          getSmallSVAssembly(singleSV, isFindLargeInsertions, assemblyData);
          return;
        }
      }
    }
  }

  // First assemble candidate contigs
  AlignData alignData;
  bool      isAssemblySuccess(false);
  isAssemblySuccess = assembleJumpContigs(
      _opt, sv, extraRefEdgeSize, extraRefSplitSize, _header, _spanningAssembler, alignData, assemblyData);
  if (!isAssemblySuccess) return;

  // Align candidate contigs back to reference
  alignJumpContigs(_opt, sv, _spanningAligner, _RNASpanningAligner, alignData, assemblyData);

  // Select the contig with the highest alignment score
  bool isContigSelected(false);
  if (_opt.isRNA) {
    isContigSelected = selectJumpContigRNA(assemblyData, _contigFilterAlignmentScores);
  } else {
    isContigSelected = selectJumpContigDNA(assemblyData, _contigFilterAlignmentScores);
  }
  if (!isContigSelected) return;

#ifdef DEBUG_REFINER
  log_os << __FUNCTION__ << ": highscoreid: " << assemblyData.bestAlignmentIndex
         << " alignment: " << assemblyData.spanningAlignments[assemblyData.bestAlignmentIndex];
#endif

  // Process the alignment into information that's easily usable
  // in the vcf output (ie. breakends in reference coordinates)
  // Summarize usable output information in a second SVBreakend
  // object -- this is the 'refined' sv:
  assemblyData.svs.push_back(sv);
  SVCandidate& newSV(assemblyData.svs.back());
  generateRefinedVCFSVCandidateFromJumpAlignment(assemblyData, newSV, _opt);

#ifdef DEBUG_REFINER
  log_os << __FUNCTION__ << ": highscore refined sv: " << newSV;
#endif
}

/// Stores information associated with each contig to assist with contig rank and selection
struct ContigScoringInfo {
  bool     isDefined   = false;
  int      score       = 0;
  unsigned index       = 0;
  unsigned variantSize = 0;
  bool     isJumped    = false;
};

void SVCandidateAssemblyRefiner::getSmallSVAssembly(
    const SVCandidate& sv, const bool isFindLargeInsertions, SVCandidateAssemblyData& assemblyData) const
{
#ifdef DEBUG_REFINER
  log_os << __FUNCTION__ << ": START\n";
#endif

  assemblyData.isSpanning = false;

  if (assemblyData.isCandidateSpanning) {
    _spanToComplexAssmRegions.addInterval(sv.bp1.interval);
  } else {
    // Check if we've already assembled this region while processing a spanning candidate, if so skip it.
    // TODO: why does this assume the spanning candidate will be run first?
    if (_spanToComplexAssmRegions.isSubsetOfRegion(sv.bp1.interval)) {
      assemblyData.isOverlapSkip = true;
      return;
    }
  }

  // how much additional reference sequence should we extract from
  // around each side of the breakend region?
  //
  // see extended description in getJumpAssembly
  static const pos_t extraRefEdgeSize(700);

  // how much reference should we additionally extract for split
  // read alignment, but not for variant-discovery alignment?
  //
  // see extended description in getJumpAssembly
  static const pos_t extraRefSplitSize(100);

  static const pos_t extraRefSize(extraRefEdgeSize + extraRefSplitSize);

  // There is a small chance that the target assembly/alignment region will not intersect with the reference.
  // This can't be treated as a bug because some chromosomes are circular. Manta is not setup to handle SVs
  // spanning the origin of circular chromosomes correctly, so these cases are skipped.
  //
  if (!isRefRegionValid(_header, sv.bp1.interval)) return;

  unsigned leadingTrim;
  unsigned trailingTrim;
  getIntervalReferenceSegment(
      _opt.referenceFilename,
      _header,
      extraRefSize,
      sv.bp1.interval,
      assemblyData.bp1ref,
      leadingTrim,
      trailingTrim);

  // The *Cut values below represent sequence which will be removed from the edges of the reference region for
  // each breakend. In most cases: (1) leadingCut and trailingCut will equal extraRefSplitSize (2)
  // maxLeadingCut and maxTrailingCut will equal extraRefSize. Sometimes these values are forced to be shorter
  // because we didn't retrieve as much reference sequence as targeted.
  const pos_t maxLeadingCut(std::max(0, extraRefSize - static_cast<pos_t>(leadingTrim)));
  const pos_t maxTrailingCut(std::max(0, extraRefSize - static_cast<pos_t>(trailingTrim)));
  const pos_t leadingCut(std::max(0, maxLeadingCut - extraRefEdgeSize));
  const pos_t trailingCut(std::max(0, maxTrailingCut - extraRefEdgeSize));

  const std::string& align1RefStr(assemblyData.bp1ref.seq());

  const bool isSearchRemoteInsertionReads(_opt.enableRemoteReadRetrieval && isFindLargeInsertions);

  // assemble contigs in the breakend region
  _smallSVAssembler.assembleComplexSVCandidate(
      sv.bp1,
      assemblyData.bp1ref,
      isSearchRemoteInsertionReads,
      assemblyData.remoteReads,
      assemblyData.contigs);

#ifdef DEBUG_REFINER
  log_os << __FUNCTION__ << ": align1RefSize/Seq: " << align1RefStr.size() << '\n';
  printSeq(align1RefStr, log_os);
  log_os << '\n';
#endif

  const unsigned contigCount(assemblyData.contigs.size());

#ifdef DEBUG_REFINER
  log_os << __FUNCTION__ << ": contigCount: " << contigCount << '\n';
  for (unsigned contigIndex(0); contigIndex < contigCount; ++contigIndex) {
    const AssembledContig& contig(assemblyData.contigs[contigIndex]);
    log_os << __FUNCTION__ << ": contigIndex: " << contigIndex << " contig: " << contig;
  }
#endif

  // make sure an alignment object exists for every contig, even if
  // it's empty:
  assemblyData.smallSVAlignments.resize(contigCount);
  assemblyData.smallSVSegments.resize(contigCount);
  assemblyData.largeInsertInfo.resize(contigCount);
  assemblyData.extendedContigs.resize(contigCount);

  // track best and second-best contig details:
  ContigScoringInfo rank1Contig;
  ContigScoringInfo rank2Contig;

  std::vector<unsigned> largeInsertionCandidateIndex;

  for (unsigned contigIndex(0); contigIndex < contigCount; ++contigIndex) {
    const AssembledContig& contig(assemblyData.contigs[contigIndex]);

#ifdef DEBUG_REFINER
    log_os << __FUNCTION__ << ": start aligning contigIndex: " << contigIndex << '\n';
#endif

    // Run a universal aligner that can find both large events (with free gap-extension after jump)
    // and smaller events (without jump)
    SVCandidateAssemblyData::SmallAlignmentResultType& alignment(assemblyData.smallSVAlignments[contigIndex]);
    std::string&                                extendedContig(assemblyData.extendedContigs[contigIndex]);
    std::vector<std::pair<unsigned, unsigned>>& candidateSegments(assemblyData.smallSVSegments[contigIndex]);
    candidateSegments.clear();

    // Remove candidate from consideration unless we find a sufficiently large indel with good flanking
    // sequence in the contig alignment
    bool isSmallSVCandidate(false);

    // To accelerate alignment, find the earliest and latest contig kmer match in the reference, then
    // limit the target reference region used for alignment to these points (but not less than the
    // original breakend region)
    //
    // New reference span reflected in adjusted*Cut values below
    pos_t adjustedLeadingCut(leadingCut);
    pos_t adjustedTrailingCut(trailingCut);
    {
      // Pick a relatively low value for k because this is just a simple runtime optimization
      static const int                merSize(10);
      std::unordered_set<std::string> contigHash;
      const unsigned                  contigSize(contig.seq.size());
      for (unsigned contigMerIndex(0); contigMerIndex < (contigSize - (merSize - 1)); ++contigMerIndex) {
        contigHash.insert(contig.seq.substr(contigMerIndex, merSize));
      }

      const pos_t refSize(align1RefStr.size());
      const pos_t minRefIndex(leadingCut);
      const pos_t maxRefIndex(refSize - (trailingCut + merSize));

      const pos_t maxFwdRefIndex(std::min(maxLeadingCut, maxRefIndex));
      pos_t       refIndex = minRefIndex;
      for (refIndex = minRefIndex; refIndex <= maxFwdRefIndex; refIndex++) {
        if (contigHash.count(align1RefStr.substr(refIndex, merSize)) != 0) break;
      }
      adjustedLeadingCut = refIndex;

      const pos_t minRevRefIndex(std::max(minRefIndex, refSize - maxTrailingCut));
      for (refIndex = (maxRefIndex); refIndex >= minRevRefIndex; refIndex--) {
        if (contigHash.count(align1RefStr.substr(refIndex, merSize)) != 0) break;
      }
      adjustedTrailingCut = (refSize - (refIndex + merSize));
    }

#ifdef DEBUG_REFINER
    log_os << __FUNCTION__ << " Ref for alignment: "
           << align1RefStr.substr(
                  adjustedLeadingCut, align1RefStr.size() - adjustedLeadingCut - adjustedTrailingCut)
           << '\n';
#endif

    // sanity check aligner requirements
    if (contig.seq.empty()) {
      using namespace illumina::common;

      std::ostringstream oss;
      oss << "Assembly produced unexpected zero-length contig. ContigIndex: " << contigIndex
          << " ContigCount: " << contigCount << " ContigDetails: " << contig;
      BOOST_THROW_EXCEPTION(GeneralException(oss.str()));
    }

    // Use a universal aligner to discover any small/large indels in the contig/reference alignment.
    {
      _largeSVAligner.align(
          contig.seq.begin(),
          contig.seq.end(),
          align1RefStr.begin() + adjustedLeadingCut,
          align1RefStr.end() - adjustedTrailingCut,
          alignment);
      alignment.align.beginPos += adjustedLeadingCut;
      getExtendedContig(alignment, contig.seq, align1RefStr, extendedContig);

      // Test alignment quality to determine if it nominates an indel candidate.
      //
      // Test is run twice, with different flanking test sizes, this way we
      // account for multiple neighboring noise scenarios.
      //
      static const unsigned spanSet[] = {100, 200};
      for (const unsigned maxQCRefSpan : spanSet) {
        std::vector<std::pair<unsigned, unsigned>> segments;
        const bool isCandidate(findCandidateVariantsFromComplexSVContigAlignment(
            maxQCRefSpan,
            _contigFilterAlignmentScores,
            alignment.align,
            contig.seq,
            align1RefStr,
            _opt.scanOpt.minCandidateVariantSize,
            segments));

        if (isCandidate) {
          // in case both ref spans are accepted take the one with the larger segment count:
          if (segments.size() > candidateSegments.size()) {
            candidateSegments = segments;
          }
          // Set isSmallSVCandidate true if an indel candidate can be nominated from alignment
          isSmallSVCandidate = true;
        }
      }

#ifdef DEBUG_REFINER
      log_os << __FUNCTION__ << ": finished aligner. contigIndex: " << contigIndex
             << " isSmallSVCandidate: " << isSmallSVCandidate << " alignment: " << alignment;
#endif
    }

    // Test each alignment for suitability to be the left or right side of a large insertion.
    //
    // All practical combinations of left and right candidates will be enumerated below to see
    // if any fit a large insertion hypothesis.
    //
    if (isFindLargeInsertions) {
      LargeInsertionInfo& candidateInsertInfo(assemblyData.largeInsertInfo[contigIndex]);
      candidateInsertInfo.clear();

      LargeInsertionInfo insertInfo;
#ifdef DEBUG_CONTIG
      log_os << __FUNCTION__ << ": contig length: " << contig.seq.size() << "\n"
             << __FUNCTION__ << ": contig seq: " << contig.seq << "\n"
             << __FUNCTION__ << ": trim contig start offset: " << contig.conservativeRange.begin_pos() << "\n"
             << __FUNCTION__ << ": trim contig end offset: " << contig.conservativeRange.end_pos() << "\n";
#endif
      ALIGNPATH::path_t apath_conservative(alignment.align.apath);
      apath_limit_read_length(contig.conservativeRange, apath_conservative);

      bool isCandidate(isLargeInsertAlignment(_largeInsertEdgeAligner, apath_conservative, insertInfo));

      if (isCandidate) {
        // if passed, then get corrected insertInfo without
        // using conservativeRange:
        LargeInsertionInfo insertInfo2;
        isCandidate = isLargeInsertAlignment(_largeInsertEdgeAligner, alignment.align.apath, insertInfo2);

        if ((insertInfo.isLeftCandidate != insertInfo2.isLeftCandidate) ||
            (insertInfo.isRightCandidate != insertInfo2.isRightCandidate)) {
          isCandidate = false;
        }

        insertInfo.contigOffset = insertInfo2.contigOffset;
        insertInfo.refOffset    = insertInfo2.refOffset;
      }

      if (isCandidate) {
        candidateInsertInfo = insertInfo;
#ifdef DEBUG_REFINER
        log_os << __FUNCTION__ << ": inserting large insertion candidation: " << candidateInsertInfo << "\n";
#endif
        largeInsertionCandidateIndex.push_back(contigIndex);
      }
    }

    if (isSmallSVCandidate) {
      /// Updates \p contigInfo to contain data for the current contig
      auto refreshContigScoringInfo = [&](ContigScoringInfo& contigInfo) {
        contigInfo.isDefined   = true;
        contigInfo.index       = contigIndex;
        contigInfo.score       = alignment.score;
        contigInfo.variantSize = getLargestIndelSize(alignment.align.apath, candidateSegments);
        contigInfo.isJumped    = alignment.isJumped;
      };

      // keep the top two highest scoring QC'd candidate:
      // TODO: we should keep all QC'd candidates for the small event case
      // FIXME : prevents us from finding overlapping events, keep vector of high-scoring contigs?

      // The current contig is selected if any of the following criteria is met:
      // (1) no contig is selected yet
      // (2) the current contig contains JUMP/JUMPINS state, but rank1Contig does not
      // (3) both (a) and (b) are true:
      //     (a) if both contigs contain JUMP/JUMPINS state (bothJumped), or
      //         if both contigs do not contain JUMP/JUMPINS state (bothNotJumped)
      //     (b) score of the current contig is higher than that of the rank1Contig
      const bool bothJumped    = alignment.isJumped && rank1Contig.isJumped;
      const bool bothNotJumped = (!alignment.isJumped) && (!rank1Contig.isJumped);

      if ((!rank1Contig.isDefined) || (alignment.isJumped && (!rank1Contig.isJumped)) ||
          (((bothJumped || bothNotJumped) && (alignment.score > rank1Contig.score)))) {
        if (rank1Contig.isDefined) {
          rank2Contig = rank1Contig;

#ifdef DEBUG_REFINER
          log_os << __FUNCTION__ << ": contigIndex: " << rank2Contig.index << " is the second high score\n";
#endif
        }

        refreshContigScoringInfo(rank1Contig);

#ifdef DEBUG_REFINER
        log_os << __FUNCTION__ << ": contigIndex: " << rank1Contig.index << " is high score\n";
#endif
      } else if ((!rank2Contig.isDefined) || (alignment.score > rank2Contig.score)) {
        refreshContigScoringInfo(rank2Contig);

#ifdef DEBUG_REFINER
        log_os << __FUNCTION__ << ": contigIndex: " << rank2Contig.index << " is the second high score\n";
#endif
      }
    }
  }

  // Test to see if the second-best contig alignment should be preferred over the best contig alignment:
  if (rank2Contig.isDefined) {
    assert(rank1Contig.isDefined);

    const unsigned rank1ContigSupportReadCount = assemblyData.contigs[rank1Contig.index].supportReads.size();
    const unsigned rank2ContigSupportReadCount = assemblyData.contigs[rank2Contig.index].supportReads.size();

#ifdef DEBUG_REFINER
    log_os << __FUNCTION__ << ": contig #" << rank1Contig.index << " has " << rank1ContigSupportReadCount
           << " support reads, with max variant size " << rank1Contig.variantSize
           << ", contains JUMP/JUMPINS state " << rank1Contig.isJumped << "\n";
    log_os << __FUNCTION__ << ": contig #" << rank2Contig.index << " has " << rank2ContigSupportReadCount
           << " support reads, with max variant size " << rank2Contig.variantSize
           << ", contains JUMP/JUMPINS state " << rank2Contig.isJumped << "\n";
#endif

    static const float minScoreRatio(0.9f);
    static const float minSupportReadCountRatio(1.2f);
    static const float minVariantSizeRatio(1.1f);

    // rank1Contig will be selected, if rank1Contig contains JUMP/JUMPINS state
    // and rank2Contig does not
    const bool rank1IsSelected((rank1Contig.isJumped && (!rank2Contig.isJumped)));
    if (!rank1IsSelected) {
      // At this point, only possibilities are:
      // both contigs contain JUMP/JUMPINS state OR
      // both contigs do not contain JUMP/JUMPINS state
      //
      // The reason is that the previous contig selection process already
      // eliminates the possibility in which
      // rank1Contig does NOT contain a jump state while rank2Contig does
      //
      // The second best contig is selected if both (1) and (2) are true:
      // (1) score2/score1 is higher than minScoreRatio
      // (2) either (a) or (b) is true:
      //     (a) supportReadCount2/supportReadCount1 is higher than minSupportReadCountRatio
      //     (b) variantSize2/variantSize1 is higher than minVariantSizeRatio
      //
      const bool rank2IsBest(
          (rank2Contig.score > (rank1Contig.score * minScoreRatio)) &&
          ((rank2ContigSupportReadCount > (rank1ContigSupportReadCount * minSupportReadCountRatio)) ||
           (rank2Contig.variantSize > (rank1Contig.variantSize * minVariantSizeRatio))));

      if (rank2IsBest) {
        rank1Contig = rank2Contig;
      }
    }

#ifdef DEBUG_REFINER
    log_os << __FUNCTION__ << ": contigIndex: " << rank1Contig.index << " is finally selected.\n";
#endif
  }

  // set any additional QC steps before deciding an alignment is usable:
  // TODO:
  std::set<pos_t> insPos;

  // finished QC, skip small deletions if no candidates have appeared:
  if (rank1Contig.isDefined) {
    assemblyData.bestAlignmentIndex = rank1Contig.index;
#ifdef DEBUG_REFINER
    log_os << __FUNCTION__ << ": highscoreid: " << rank1Contig.index
           << " alignment: " << assemblyData.smallSVAlignments[rank1Contig.index];
#endif

    // process the alignment into information that's easily usable
    // in the vcf output (ie. breakends in reference coordinates)

    const AssembledContig& bestContig(assemblyData.contigs[assemblyData.bestAlignmentIndex]);
    const SVCandidateAssemblyData::SmallAlignmentResultType& bestAlign(
        assemblyData.smallSVAlignments[assemblyData.bestAlignmentIndex]);

    const SVCandidateAssemblyData::CandidateSegmentSetType& candidateSegments(
        assemblyData.smallSVSegments[assemblyData.bestAlignmentIndex]);
    unsigned segmentIndex = 0;
    for (const SVCandidateAssemblyData::CandidateSegmentType& segRange : candidateSegments) {
      // copy the low-res candidate sv and start customizing:
      assemblyData.svs.push_back(sv);

      SVCandidate& newSV(assemblyData.svs.back());
      newSV.assemblyAlignIndex   = assemblyData.bestAlignmentIndex;
      newSV.assemblySegmentIndex = segmentIndex;
      setSmallCandSV(assemblyData.bp1ref, bestContig.seq, bestAlign.align, segRange, newSV, _opt);
      segmentIndex++;

      // provide a weak filter to keep fully and partially
      // assembled duplicates of the same event from occurring:
      if (getExtendedSVType(newSV) == EXTENDED_SV_TYPE::INSERT) {
        insPos.insert(newSV.bp1.interval.range.begin_pos());
      }

#ifdef DEBUG_REFINER
      log_os << __FUNCTION__ << ": small refined sv: " << newSV;
#endif

#ifdef DEBUG_CONTIG
      const int               contigSize = bestContig.seq.length();
      const ALIGNPATH::path_t apathTillSvStart(
          &bestAlign.align.apath[0], &bestAlign.align.apath[segRange.first]);
      const ALIGNPATH::path_t apathTillSvEnd(
          &bestAlign.align.apath[0], &bestAlign.align.apath[segRange.second + 1]);
      const int leftSize  = apath_read_length(apathTillSvStart);
      const int endPos    = apath_read_length(apathTillSvEnd);
      const int rightSize = contigSize - apath_read_length(apathTillSvEnd);

      log_os << __FUNCTION__ << ": contig has size " << contigSize << ": " << bestContig.seq << "\n";
      log_os << __FUNCTION__ << ": left part has size " << leftSize << ": "
             << bestContig.seq.substr(0, leftSize) << "\n";
      log_os << __FUNCTION__ << ": right part has size " << rightSize << ": "
             << bestContig.seq.substr(endPos, rightSize) << "\n";
#endif
    }
  }

  // search for large deletions with incomplete insertion sequence
  // assembly:
  {
    // In case of no fully-assembled candidate, solve for any
    // strong large insertion candidate
    //
    // This is done by searching through combinations of the left
    // and right insertion side candidates found in the primary
    // contig processing loop
    if (isFindLargeInsertions) {
      processLargeInsertion(
          sv,
          leadingCut,
          trailingCut,
          _largeInsertCompleteAligner,
          largeInsertionCandidateIndex,
          insPos,
          assemblyData,
          _opt);
    }
  }
}
