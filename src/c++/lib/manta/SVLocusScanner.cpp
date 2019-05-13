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
/// \author Chris Saunders
/// \author Ole Schulz-Trieglaff
/// \author Bret Barnes
///

#include "manta/SVLocusScanner.hpp"
#include "blt_util/align_path_util.hpp"
#include "blt_util/parse_util.hpp"
#include "blt_util/string_util.hpp"
#include "common/Exceptions.hpp"
#include "htsapi/SimpleAlignment_bam_util.hpp"
#include "htsapi/align_path_bam_util.hpp"
#include "htsapi/bam_record_util.hpp"
#include "manta/RemoteMateReadUtil.hpp"
#include "manta/SVCandidateUtil.hpp"
#include "manta/SVLocusScannerSemiAligned.hpp"

#include <iostream>

//#define DEBUG_SCANNER

//#define DEBUG_IS_SHADOW

#ifdef DEBUG_SCANNER
#include "blt_util/log.hpp"
#endif

/// \brief Utilities pertaining to classifying anomolous fragments based on size so that they can be
/// treated/weighted differently.
///
namespace FragmentSizeType {
/// \brief Fragments within this factor of the minimum size cutoff are treated as 'close' pairs and receive a
/// modified evidence count.
static const float closePairFactor(4);

static const float minLargeEventRegionFactor(10);
static const float maxLargeEventRegionFactor(20);

static index_t classifySize(const SVLocusScanner::CachedReadGroupStats& rgStats, const int fragmentSize)
{
  if (fragmentSize < rgStats.properPair.min) return COMPRESSED;
  if (fragmentSize > rgStats.properPair.max) {
    if (fragmentSize < rgStats.minDistantFragmentSize) return CLOSE;
    return DISTANT;
  }
  return NORMAL;
}

static bool isLarge(const index_t i)
{
  switch (i) {
  case NORMAL:
  case COMPRESSED:
    return false;
  default:
    return true;
  }
}
}  // namespace FragmentSizeType

/// \brief Convert input details into an SVObservation indicating an SV candidate which occurs on a single
/// chromosome from \p leftPos to \p rightPos.
///
/// \param[in] candidateChromosomeIndex Index of the chromosome on which the SV candidate occurs
///
/// \param[in] svEvidenceSource The type of SV evidence (anomalous read pair, CIGAR string, etc)
///
/// \param[in] dnaFragmentSVEvidenceSource The source of SV evidence within the read fragment (read1, read2,
/// etc..)
///
/// \param[in] isComplex If true, create a 'complex' candidate indicating a non-specific local assembly task
/// to search for indels in a single region.
static SVObservation getSplitSVCandidate(
    const ReadScannerDerivOptions&                  dopt,
    const int32_t                                   candidateChromosomeIndex,
    const pos_t                                     candidateChromosomeLength,
    const pos_t                                     leftPos,
    const pos_t                                     rightPos,
    const SVEvidenceType::index_t&                  svEvidenceSource,
    const SourceOfSVEvidenceInDNAFragment::index_t& dnaFragmentSVEvidenceSource,
    const bool                                      isComplex = false)
{
  SVObservation sv;
  SVBreakend&   localBreakend(sv.bp1);
  SVBreakend&   remoteBreakend(sv.bp2);

  localBreakend.interval.tid  = candidateChromosomeIndex;
  remoteBreakend.interval.tid = candidateChromosomeIndex;

  localBreakend.lowresEvidence.add(svEvidenceSource);
  sv.svEvidenceType              = svEvidenceSource;
  sv.dnaFragmentSVEvidenceSource = dnaFragmentSVEvidenceSource;

  if (!isComplex) {
    remoteBreakend.lowresEvidence.add(svEvidenceSource);
    localBreakend.state  = SVBreakendState::RIGHT_OPEN;
    remoteBreakend.state = SVBreakendState::LEFT_OPEN;
  } else {
    localBreakend.state  = SVBreakendState::COMPLEX;
    remoteBreakend.state = SVBreakendState::UNKNOWN;
  }

  localBreakend.interval.range.set_begin_pos(std::max(0, leftPos - dopt.beforeBreakend));

  if (!isComplex) {
    localBreakend.interval.range.set_end_pos(
        std::min(candidateChromosomeLength, (leftPos + dopt.afterBreakend)));
  } else {
    localBreakend.interval.range.set_end_pos(
        std::min(candidateChromosomeLength, (rightPos + dopt.afterBreakend)));
  }

  remoteBreakend.interval.range.set_begin_pos(std::max(0, rightPos - dopt.beforeBreakend));
  remoteBreakend.interval.range.set_end_pos(
      std::min(candidateChromosomeLength, (rightPos + dopt.afterBreakend)));

  return sv;
}

/// determine, based on clipping in the cigar string, if this split alignment
/// has its breakpoint on the downstream (right) end or the upstream (left) end
static bool isSplitOpenDownstream(const ALIGNPATH::path_t& align)
{
  using namespace ALIGNPATH;
  /// TODO replace this heuristic with a better check (looking at all SA alignments at once)
  return (apath_clip_lead_size(align) < apath_clip_trail_size(align));
}

static void updateSABreakend(
    const ReadScannerDerivOptions& dopt,
    const SimpleAlignment&         align,
    SVBreakend&                    breakend,
    const bam_header_info&         bamHeader)
{
  // Need to use the match descriptors to determine if the split is upstream (i.e. 5' assuming fwd strand)
  // of the current alignment (i.e. we are clipped on the left side) or downstream
  // Below is the logic to convert these  to breakend candidates (everything is relative to the forward
  // strand):
  //
  // DownStream => RIGHT_OPEN
  // Upstream => LEFT_OPEN
  //

  const bool isSplitDownstream(isSplitOpenDownstream(align.path));

  if (isSplitDownstream) {
    breakend.state = SVBreakendState::RIGHT_OPEN;
  } else {
    breakend.state = SVBreakendState::LEFT_OPEN;
  }

  breakend.interval.tid = align.tid;
  const pos_t alignTlength(bamHeader.chrom_data[align.tid].length);
  // get the position of the breakend implied by the split, if split
  // is downstream (see above) the split position is the end of this split
  // read segment
  int pos = align.pos;
  if (isSplitDownstream) {
    using namespace ALIGNPATH;
    pos += apath_ref_length(align.path);
  }
  breakend.interval.range.set_begin_pos(std::max(0, pos - dopt.beforeBreakend));
  breakend.interval.range.set_end_pos(std::min(alignTlength, (pos + dopt.afterBreakend)));
}

/// \brief Convert details from one of \p localRead's SA-tag split-read alignments into an SVObservation.
///
/// \param[in] dnaFragmentSVEvidenceSource The source of SV evidence within the read fragment (read1, read2,
/// etc..)
static SVObservation getSplitSACandidate(
    const ReadScannerDerivOptions&                 dopt,
    const bam_record&                              localRead,
    const SimpleAlignment&                         localAlign,
    const SimpleAlignment&                         remoteAlign,
    const SourceOfSVEvidenceInDNAFragment::index_t dnaFragmentSVEvidenceSource,
    const bam_header_info&                         bamHeader)
{
  using namespace SVEvidenceType;
  static const index_t svSource(SPLIT_ALIGN);

  SVObservation sv;
  sv.svEvidenceType              = svSource;
  sv.dnaFragmentSVEvidenceSource = dnaFragmentSVEvidenceSource;

  SVBreakend& localBreakend(sv.bp1);
  SVBreakend& remoteBreakend(sv.bp2);

  // use single-side evidence, have to read the supp read to get the
  // reverse edge. this protects against double-count:
  localBreakend.lowresEvidence.add(svSource);

  updateSABreakend(dopt, localAlign, localBreakend, bamHeader);
  updateSABreakend(dopt, remoteAlign, remoteBreakend, bamHeader);

  // If the local (bp1) alignment is split downstream (on the right side) then this read goes from bp1 -> bp2.
  // If it is a forward read (e.g. read1 on + strand), this means it's a forward read for this event.
  const bool isSplitDownstream(isSplitOpenDownstream(localAlign.path));
  const bool isReadFw = (localRead.is_first() == localRead.is_fwd_strand());
  if (dopt.isTranscriptStrandKnown) {
    if (isReadFw == isSplitDownstream) {
      sv.forwardTranscriptStrandReadCount += 1;
    } else {
      sv.reverseTranscriptStrandReadCount += 1;
    }
  }
  return sv;
}

typedef std::map<std::string, int32_t> chromMap_t;

/// \brief Enumerate the fields for each split alignment segment as described in the BAM spec "SA" tag entry
namespace SAFields {
enum index_t { CHROM, POS, STRAND, CIGAR, MAPQ, NM, SIZE };
}

static void parseSACandidatesFromRead(
    const ReadScannerOptions&     opt,
    const bam_record&             bamRead,
    const chromMap_t&             chromToIndex,
    std::vector<SimpleAlignment>& splitAlignments)
{
  using namespace illumina::common;
  using namespace ALIGNPATH;

  splitAlignments.clear();

  std::vector<std::string> splitAlignmentSegmentStrings;
  {
    static const char splitAlignmentTag[] = {'S', 'A'};
    const char*       splitAlignmentString(bamRead.get_string_tag(splitAlignmentTag));
    if (nullptr == splitAlignmentString) return;

    split_string(splitAlignmentString, ';', splitAlignmentSegmentStrings);
    if ((!splitAlignmentSegmentStrings.empty()) && splitAlignmentSegmentStrings.back().empty()) {
      splitAlignmentSegmentStrings.pop_back();
    }
  }

  // Only handle a single split alignment segment right now (meaning the read is split into 2 parts).
  // In the future we could sort the SA tags by order on the template, possibly
  // also removing segments that map to two different areas.
  if (splitAlignmentSegmentStrings.size() > 1) return;

  for (const std::string& splitAlignmentSegmentString : splitAlignmentSegmentStrings) {
#ifdef DEBUG_SCANNER
    log_os << __FUNCTION__ << ": splitAlignmentSegmentString: " << splitAlignmentSegmentString << "\n";
#endif
    std::vector<std::string> splitAlignmentSegmentFields;
    split_string(splitAlignmentSegmentString, ',', splitAlignmentSegmentFields);

    if (splitAlignmentSegmentFields.size() != SAFields::SIZE) {
      std::ostringstream oss;
      oss << "Unexpected format in the split alignment segment: '" << splitAlignmentSegmentString << "'";
      BOOST_THROW_EXCEPTION(GeneralException(oss.str()));
    }

    /// filter split reads with low MappingQuality:
    const unsigned splitAlignmentMappingQuality(
        illumina::blt_util::parse_unsigned_str(splitAlignmentSegmentFields[SAFields::MAPQ]));
    if (splitAlignmentMappingQuality < opt.minMapq) continue;

    const std::string&               splitAlignmentChrom(splitAlignmentSegmentFields[SAFields::CHROM]);
    const chromMap_t::const_iterator ci(chromToIndex.find(splitAlignmentChrom));

    if (ci == chromToIndex.end()) {
      std::ostringstream oss;
      oss << "Split alignment segment maps to an unknown chromosome: '" << splitAlignmentChrom << "'";
      BOOST_THROW_EXCEPTION(GeneralException(oss.str()));
    }

    splitAlignments.emplace_back();
    SimpleAlignment& sal(splitAlignments.back());
    sal.tid = (ci->second);  // convert chr to int32_t via new bam header map
    sal.pos = (illumina::blt_util::parse_int_str(splitAlignmentSegmentFields[SAFields::POS]) - 1);
    {
      const char splitAlignmentStrand(splitAlignmentSegmentFields[SAFields::STRAND][0]);  // convert to char
      if (!((splitAlignmentStrand == '-') || (splitAlignmentStrand == '+'))) {
        std::ostringstream oss;
        oss << "Unexpected strand entry in split alignment segment: '" << splitAlignmentSegmentString << "'";
        BOOST_THROW_EXCEPTION(GeneralException(oss.str()));
      }
      sal.is_fwd_strand = (splitAlignmentStrand == '+');
    }

    cigar_to_apath(splitAlignmentSegmentFields[SAFields::CIGAR].c_str(), sal.path);
  }
}

/// \brief Find all 'SA' formatted split read sub-alignments in a single read alignment record, and convert
/// these into SVObservation objects.
///
/// Convert each BAM "SA" sub-alignment record in a read into an SVObservation object. The method currently
/// only returns an SVObservation object if the read contains exactly one "SA" sub-alignment record, but the
/// interface supports parsing out any number of split alignments, and this function may be updated to do so
/// in the future.
///
/// \param[in] dnaFragmentSVEvidenceSource The source of SV evidence within the read fragment (read1, read2,
/// etc..)
///
/// \param[in,out] candidates New SVObservation objects are appended to this vector. Contents of the vector
/// are preserved but not read.
static void getSACandidatesFromRead(
    const ReadScannerOptions&                      opt,
    const ReadScannerDerivOptions&                 dopt,
    const bam_record&                              localRead,
    const SimpleAlignment&                         localAlign,
    const SourceOfSVEvidenceInDNAFragment::index_t dnaFragmentSVEvidenceSource,
    const bam_header_info&                         bamHeader,
    std::vector<SVObservation>&                    candidates)
{
  using namespace ALIGNPATH;

  std::vector<SimpleAlignment> remoteAlign;
  parseSACandidatesFromRead(opt, localRead, bamHeader.chrom_to_index, remoteAlign);

  if (remoteAlign.empty()) return;

  // Only handle a single split alignment right now. In the future we should sort the SA tags by order on the
  // template, possibly also removing segments that map to two different areas.
  if (remoteAlign.size() > 1) return;

  for (const auto& ral : remoteAlign) {
    candidates.push_back(
        getSplitSACandidate(dopt, localRead, localAlign, ral, dnaFragmentSVEvidenceSource, bamHeader));
#ifdef DEBUG_SCANNER
    log_os << __FUNCTION__ << ": evaluating SA sv for inclusion: " << candidates.back() << "\n";
#endif
  }
}

/// \brief Convert all large indels already present in a single read's alignment into SVObservation objects
///
/// \param[in] dnaFragmentSVEvidenceSource The source of SV evidence within the read fragment (read1, read2,
/// etc..) \param[in,out] candidates New SVObservation objects are appended to this vector. Contents of the
/// vector are preserved but not read.
static void getSVCandidatesFromReadIndels(
    const ReadScannerOptions&                      opt,
    const ReadScannerDerivOptions&                 dopt,
    const SimpleAlignment&                         align,
    const bam_header_info&                         bamHeader,
    const SourceOfSVEvidenceInDNAFragment::index_t dnaFragmentSVEvidenceSource,
    std::vector<SVObservation>&                    candidates)
{
  using namespace SVEvidenceType;
  static const index_t svSource(CIGAR);

  using namespace ALIGNPATH;
  const std::pair<unsigned, unsigned> ends(get_match_edge_segments(align.path));

  unsigned      pathIndex(0);
  unsigned      readOffset(0);
  pos_t         refHeadPos(align.pos);
  const int32_t chromosomeIndex(align.tid);
  const pos_t   chromosomeLength(bamHeader.chrom_data[chromosomeIndex].length);

  const unsigned pathSize(align.path.size());
  while (pathIndex < pathSize) {
    const path_segment& ps(align.path[pathIndex]);
    const bool          isBeginEdge(pathIndex < ends.first);
    const bool          isEndEdge(pathIndex > ends.second);
    const bool          isEdgeSegment(isBeginEdge || isEndEdge);

    // in this case, swap means combined insertion/deletion
    const bool isSwapStart(is_segment_swap_start(align.path, pathIndex));

    if (isEdgeSegment && isSwapStart) {
      using namespace illumina::common;

      std::ostringstream oss;
      oss << "Can't process unexpected alignment pattern: " << align << "\n";
      BOOST_THROW_EXCEPTION(GeneralException(oss.str()));
    }

    unsigned nPathSegments(1);  // number of path segments consumed
    if (isEdgeSegment) {
      // edge inserts are allowed for intron adjacent and grouper reads, edge deletions for intron adjacent
      // only

      if (ps.type == INSERT) {
        if (ps.length >= opt.minCandidateVariantSize) {
          static const bool isComplex(true);
          candidates.push_back(getSplitSVCandidate(
              dopt,
              chromosomeIndex,
              chromosomeLength,
              refHeadPos,
              refHeadPos,
              svSource,
              dnaFragmentSVEvidenceSource,
              isComplex));
        }
      }
    } else if (isSwapStart) {
      const swap_info sinfo(align.path, pathIndex);
      if ((sinfo.delete_length >= opt.minCandidateVariantSize) ||
          (sinfo.insert_length >= opt.minCandidateVariantSize)) {
        candidates.push_back(getSplitSVCandidate(
            dopt,
            chromosomeIndex,
            chromosomeLength,
            refHeadPos,
            refHeadPos + sinfo.delete_length,
            svSource,
            dnaFragmentSVEvidenceSource));
      }

      nPathSegments = sinfo.n_seg;
    } else if (is_segment_type_indel(align.path[pathIndex].type)) {
      // regular indel:

      if (ps.type == DELETE) {
        if (ps.length >= opt.minCandidateVariantSize) {
          candidates.push_back(getSplitSVCandidate(
              dopt,
              chromosomeIndex,
              chromosomeLength,
              refHeadPos,
              refHeadPos + ps.length,
              svSource,
              dnaFragmentSVEvidenceSource));
        }
      } else if (ps.type == INSERT) {
        if (ps.length >= opt.minCandidateVariantSize) {
          candidates.push_back(getSplitSVCandidate(
              dopt,
              chromosomeIndex,
              chromosomeLength,
              refHeadPos,
              refHeadPos,
              svSource,
              dnaFragmentSVEvidenceSource));
        }
      }
    }

    for (unsigned i(0); i < nPathSegments; ++i) {
      increment_path(align.path, pathIndex, readOffset, refHeadPos);
    }
  }
}

/// \brief Detect if \p bamRead is semi-aligned (has one or two poorly-aligned ends), if so convert each
/// poorly aligned read end to an SVObservation object.
///
/// Reads with poorly aligned ends are referred to as semi-aligned because part of the read is confidently
/// aligned and part of the read (one or both edges) is poorly aligned. The basis for an SV hypothesis is that
/// the poorly aligned sections of the reads should be mapped elsewhere due to a large indel or SV breakpoint
/// not reflected in the current read alignment.
///
/// \param[in] dnaFragmentSVEvidenceSource The source of SV evidence within the read fragment (read1, read2,
/// etc..)
///
/// \param[in,out] candidates New SVObservation objects are appended to this vector. Contents of the vector
/// are preserved but not read.
static void getSVCandidatesFromSemiAligned(
    const ReadScannerOptions&                      opt,
    const ReadScannerDerivOptions&                 dopt,
    const bam_record&                              bamRead,
    const bam_header_info&                         bamHeader,
    const SimpleAlignment&                         bamAlign,
    const SourceOfSVEvidenceInDNAFragment::index_t dnaFragmentSVEvidenceSource,
    const reference_contig_segment&                refSeq,
    std::vector<SVObservation>&                    candidates)
{
  unsigned leadingMismatchLen(0);
  unsigned trailingMismatchLen(0);
  pos_t    leadingRefPos(0), trailingRefPos(0);
  getSVBreakendCandidateSemiAligned(
      bamRead,
      bamAlign,
      refSeq,
      opt.useOverlapPairEvidence,
      leadingMismatchLen,
      leadingRefPos,
      trailingMismatchLen,
      trailingRefPos);

  if ((leadingMismatchLen + trailingMismatchLen) >= bamRead.read_size()) return;

  using namespace SVEvidenceType;
  static const index_t svEvidenceSource(SEMIALIGN);

  // semi-aligned reads don't define a full hypothesis, so they're always evidence for a 'complex' ie.
  // undefined, event in a fashion analogous to clipped reads
  static const bool isComplex(true);
  const int32_t     candidateChromosomeIndex(bamRead.target_id());
  const pos_t       candidateChromosomeLength(bamHeader.chrom_data[candidateChromosomeIndex].length);

  if (leadingMismatchLen >= opt.minSemiAlignedMismatchLen) {
    const pos_t pos(leadingRefPos);
    candidates.push_back(getSplitSVCandidate(
        dopt,
        candidateChromosomeIndex,
        candidateChromosomeLength,
        pos,
        pos,
        svEvidenceSource,
        dnaFragmentSVEvidenceSource,
        isComplex));
  }

  if (trailingMismatchLen >= opt.minSemiAlignedMismatchLen) {
    const pos_t pos(trailingRefPos);
    candidates.push_back(getSplitSVCandidate(
        dopt,
        candidateChromosomeIndex,
        candidateChromosomeLength,
        pos,
        pos,
        svEvidenceSource,
        dnaFragmentSVEvidenceSource,
        isComplex));
  }
}

namespace {

/// \brief Utility class for ::getSVCandidatesFromPair to evaluate the anomalous status of a read pair and
/// generate the corresponding SVCandidate object for anomalous cases.
///
/// The point of moving this logic to a utility class is to improve efficiency by reusing some intermediates
/// between the tests and SVCandidate generation.
///
/// Usage:
/// The class must be reset() with read pair alignment information before any other methods can be called.
///
struct AlignmentPairAnalyzer {
  AlignmentPairAnalyzer(
      const ReadScannerOptions&                   opt,
      const ReadScannerDerivOptions&              dopt,
      const SVLocusScanner::CachedReadGroupStats& rstats,
      const bam_header_info&                      bamHeader)
    : _opt(opt), _dopt(dopt), _rstats(rstats), _header(bamHeader)
  {
  }

  /// \brief Reset read pair (local and remote) alignment information
  ///
  /// This method must be called before any other methods.
  ///
  /// \param[in] localAlignment Alignment of the local read from a read pair
  ///
  /// \param[in] remoteAlignment Alignment of the remote read from a read pair
  ///
  /// \param[in] isRemoteAlignmentInferred True if remote read alignment was inferred from the local read
  /// alignment instead of directly observed from the remote read alignment record
  ///
  /// \param[in] isForwardStrand True if the local read is 1st in read pair (this value is used in stranded
  /// RNA mode)
  void reset(
      const SimpleAlignment& localAlignment,
      const SimpleAlignment& remoteAlignment,
      const bool             isRemoteAlignmentInferred,
      const bool             isForwardStrand)
  {
    _localAlignment            = &localAlignment;
    _remoteAlignment           = &remoteAlignment;
    _isRemoteAlignmentInferred = isRemoteAlignmentInferred;
    _isForwardStrand           = isForwardStrand;
    _isBreakendRegionScaleSet  = false;
    _breakendRegionScale       = 0;
    _totalNonInsertSize        = 0;
    _localEndRefPos            = 0;
    _remoteEndRefPos           = 0;
  }

  /// \brief Test if read pair is anomalous
  ///
  /// Requires that ::reset has been called at least once on this class to initialize the local alignment.
  ///
  /// An important side effect of this test is that the read pair's breakend scaling is computed in order to
  /// respond to this query. If the read is anomalous, this breakend scaling information can be efficiently
  /// reused to generate the intended breakend scaling when the pair is converted to an SV candidate.
  ///
  /// Breakends scales are used to reduce the breakend region size for small deletions as a noise reduction
  /// step.
  ///
  /// \return True if read pair is anomalous
  bool isAnomalousReadPair()
  {
    assert(isInit());
    if (!_isBreakendRegionScaleSet) setLargeEventRegionScale();
    return (_breakendRegionScale >= 0.);
  }

  /// \brief Convert the input alignments (provided in ::reset) into a structural variant observation
  ///
  /// The regions where each of the two breakends are likely to be found in the event that an SV exists are
  /// computed from the two read alignments and the fragment size distribution.
  ///
  /// Requires that ::isAnomalousReadPair has already been called since the last call to ::reset.
  ///
  /// \param[out] sv The SVObservation inferred from the anomalous read pair
  void getSVObservation(SVObservation& sv)
  {
    assert(_isBreakendRegionScaleSet);
    assert((_breakendRegionScale >= 0.) && (_breakendRegionScale <= 1.));

    using namespace SVEvidenceType;
    static const index_t svLocalPair(LOCAL_PAIR);
    static const index_t svPair(PAIR);

    sv.svEvidenceType              = svLocalPair;
    sv.dnaFragmentSVEvidenceSource = SourceOfSVEvidenceInDNAFragment::READ_PAIR;

    SVBreakend& localBreakend(sv.bp1);
    SVBreakend& remoteBreakend(sv.bp2);

    localBreakend.lowresEvidence.add(svLocalPair);

    if (_dopt.isTranscriptStrandKnown) {
      if (_isForwardStrand) {
        sv.forwardTranscriptStrandReadCount++;
      } else {
        sv.reverseTranscriptStrandReadCount++;
      }
    }

    if (!_isRemoteAlignmentInferred) {
      remoteBreakend.lowresEvidence.add(svLocalPair);
      localBreakend.lowresEvidence.add(svPair);
      remoteBreakend.lowresEvidence.add(svPair);
      sv.svEvidenceType = svPair;
    }

    // Determine the maximum (implicitly non-anomalous) fragment size to use when computing breakend regions,
    // according to the _breakendRegionScale value found earlier.
    const double maxFragmentSize(
        (_breakendRegionScale * _rstats.largeScaleEventBreakendRegion.max) +
        ((1. - _breakendRegionScale) * _rstats.breakendRegion.max));

    // The size used for the breakend is the maximum fragment size minus the length of read1 and read2 we can
    // already account for. If the true fragment size is within maximum fragment size, then the breakend
    // should be found within this remaining size from the current inside edge of read1 or read2.
    const pos_t breakendSize(std::max(
        static_cast<pos_t>(_opt.minPairBreakendSize),
        static_cast<pos_t>(maxFragmentSize - _totalNonInsertSize)));

    const pos_t localStartRefPos(localAlign().pos);
    const pos_t remoteStartRefPos(remoteAlign().pos);

    localBreakend.interval.tid = localAlign().tid;
    const pos_t localChromLength(_header.chrom_data[localBreakend.interval.tid].length);

    // Expected breakpoint range is from the edge of the read alignment to the last possible breakpoint
    // position assuming the fragment size is maxFragmentSize.
    if (localAlign().is_fwd_strand) {
      localBreakend.state = SVBreakendState::RIGHT_OPEN;
      localBreakend.interval.range.set_begin_pos(std::min(localChromLength, _localEndRefPos));
      localBreakend.interval.range.set_end_pos(std::min(localChromLength, (_localEndRefPos + breakendSize)));
    } else {
      localBreakend.state = SVBreakendState::LEFT_OPEN;
      localBreakend.interval.range.set_end_pos(localStartRefPos);
      localBreakend.interval.range.set_begin_pos(std::max(0, (localStartRefPos - breakendSize)));
    }

    remoteBreakend.interval.tid = remoteAlign().tid;
    const pos_t remoteChromLength(_header.chrom_data[remoteBreakend.interval.tid].length);
    if (remoteAlign().is_fwd_strand) {
      remoteBreakend.state = SVBreakendState::RIGHT_OPEN;
      remoteBreakend.interval.range.set_begin_pos(std::min(remoteChromLength, _remoteEndRefPos));
      remoteBreakend.interval.range.set_end_pos(
          std::min(remoteChromLength, (_remoteEndRefPos + breakendSize)));
    } else {
      remoteBreakend.state = SVBreakendState::LEFT_OPEN;
      remoteBreakend.interval.range.set_end_pos(remoteStartRefPos);
      remoteBreakend.interval.range.set_begin_pos(std::max(0, (remoteStartRefPos - breakendSize)));
    }
  }

  /// \brief Test if the 'inside' end of either read in a read pair touches or goes past either edge of the
  /// chromosome it has been mapped to
  ///
  /// Given chromosome and read pairs:
  ///         |----------chrom---------------|
  /// pairA:  <--read2--|    |--read1--->
  /// pairB:                         |-read1->   <--read2--|
  /// pairC:                   |-read1->         <--read2--|
  ///
  /// Here pairs A and B, but not C, meet the requirements to be marked true by this test.
  ///
  /// Requires that ::isAnomalousReadPair has already been called since the last call to ::reset.
  ///
  /// \return True if the inside end of either read is aligned to/past the end of the chromosome
  bool isAlignedToChromEnds() const
  {
    assert(_isBreakendRegionScaleSet);

    if (localAlign().is_fwd_strand) {
      const int32_t localChromIndex(localAlign().tid);
      const pos_t   localChromLength(_header.chrom_data[localChromIndex].length);
      if (_localEndRefPos >= localChromLength) return true;
    } else {
      if (localAlign().pos <= 0) return true;
    }

    if (remoteAlign().is_fwd_strand) {
      const int32_t remoteChromIndex(remoteAlign().tid);
      const pos_t   remoteChromLength(_header.chrom_data[remoteChromIndex].length);
      if (_remoteEndRefPos >= remoteChromLength) return true;
    } else {
      if (remoteAlign().pos <= 0) return true;
    }

    return false;
  }

private:
  /// \return The length of unaligned sequence from one read on the 'inside' of the read pair
  ///
  /// \example Returned value for each read from an example read pair
  ///
  /// For read1 and read2, with soft clip segments indicated by an S:
  ///
  ///  S\----read1----->SSSS     SS<-----read2-----|SSSSS
  ///
  /// The unaligned sequence length on the 'inside' of the pair for read1 is 4
  /// The unaligned sequence length on the 'inside' of the pair for read2 is 2
  static unsigned unalignedSequenceLengthInsideReadPair(const SimpleAlignment& al)
  {
    if (al.is_fwd_strand)
      return unalignedSuffixSize(al.path);
    else
      return unalignedPrefixSize(al.path);
  }

  /// \return Size of the read, excluding any unaligned segment on the inside of the read's pair
  static unsigned getNonInsertReadSize(const SimpleAlignment& al)
  {
    const unsigned readSize(apath_read_length(al.path));
    return readSize - unalignedSequenceLengthInsideReadPair(al);
  }

  /// \return The 0-indexed reference position one base after the final mapped base of the input alignment
  static pos_t getEndRefPos(const SimpleAlignment& al) { return (al.pos + apath_ref_length(al.path)); }

  /// \brief Determine the appropriate scaling factor to be used for breakend regions of the given read pair.
  ///
  /// Different breakend sizes are used for long-range vs short-range deletions (this is a noise reduction
  /// step). Non-deletions use the default (long-range) scaling factor. This function computes the required
  /// scaling factor.
  ///
  /// The scaling factor will be set to 0 if the read pair is found to be non-anomalous, so this computation
  /// also embeds an anomalous read pair test.
  void setLargeEventRegionScale()
  {
    // _breakendRegionScale ramps from 0 to 1 as we go from short to long deletions sizes, and remains set to
    // 1 for non-deletions.
    _isBreakendRegionScaleSet = true;
    _breakendRegionScale      = 1.0;

    // Find the read size excluding soft-clip/edge-insert on the 'inside' of the read pair's fragment
    const unsigned localNonInsertSize(getNonInsertReadSize(localAlign()));
    const unsigned remoteNonInsertSize(getNonInsertReadSize(remoteAlign()));

    // Total the 'used' read span of read1 and read2 (ie. the elements of the
    // fragment that are not part of the insert between the reads)
    _totalNonInsertSize = (localNonInsertSize + remoteNonInsertSize);

    const pos_t localStartRefPos(localAlign().pos);
    const pos_t remoteStartRefPos(remoteAlign().pos);
    _localEndRefPos  = getEndRefPos(localAlign());
    _remoteEndRefPos = getEndRefPos(remoteAlign());

    // Check if fragment size is still anomalous after accounting for read alignment patterns

    // First, jump out if the reads are aligned to different chromosomes or orientations
    // (this isn't done earlier because the values above are used to generate the SVObservation
    if ((localAlign().tid != remoteAlign().tid) ||
        (localAlign().is_fwd_strand == remoteAlign().is_fwd_strand))
      return;

    // Get the range spanned by the read pair insert
    known_pos_range2 insertRange;
    if (localAlign().is_fwd_strand) {
      insertRange.set_range(_localEndRefPos, remoteStartRefPos);
    } else {
      insertRange.set_range(_remoteEndRefPos, localStartRefPos);
    }

    // Get length of fragment after accounting for any variants described directly in either read alignment.
    // Note the size of insertRange can be negative, so don't call insertRange.size()
    const pos_t cigarAdjustedFragmentSize(
        _totalNonInsertSize + (insertRange.end_pos() - insertRange.begin_pos()));

    // Identify read pairs with an 'outtie' orientation and skip these cases (so they will be counted as
    // anomalous) These are typically small fragments from FFPE samples.
    const bool isOuttieReadPair(cigarAdjustedFragmentSize < 0);

    if (isOuttieReadPair) return;

    const bool isLargeFragment(
        cigarAdjustedFragmentSize > (_rstats.properPair.max + _opt.minCandidateVariantSize));

    if (isLargeFragment) {
      _breakendRegionScale = _rstats.largeEventRegionScaler.getScale(cigarAdjustedFragmentSize);
    } else {
      _breakendRegionScale = -1.;
    }
  }

  bool isInit() const { return (_localAlignment != nullptr); }

  const SimpleAlignment& localAlign() const { return *_localAlignment; }

  const SimpleAlignment& remoteAlign() const { return *_remoteAlignment; }

  const ReadScannerOptions&                   _opt;
  const ReadScannerDerivOptions&              _dopt;
  const SVLocusScanner::CachedReadGroupStats& _rstats;
  const bam_header_info&                      _header;
  const SimpleAlignment*                      _localAlignment  = nullptr;
  const SimpleAlignment*                      _remoteAlignment = nullptr;

  /// \brief This is true if the remote alignment is inferred indirectly from the local read's alignment
  /// record
  bool _isRemoteAlignmentInferred = true;

  /// \brief True if the local read is 1st in read pair (this value is used in stranded RNA mode)
  bool _isForwardStrand = false;

  /// \brief True if _breakendRegionScale is initialized
  bool _isBreakendRegionScaleSet = false;

  /// \brief Scaling factor applied to the breakend regions of the SV candidate generated from the given read
  /// pair
  double _breakendRegionScale = 0.;

  /// \brief The total length of local and remote reads, excluding any unaligned read segments on the inside
  /// of the read pair
  unsigned _totalNonInsertSize = 0;

  /// \brief The zero-indexed reference position one base after the local read's last aligned base
  pos_t _localEndRefPos = 0;

  /// \brief The zero-indexed reference position one base after the remote read's last aligned base
  pos_t _remoteEndRefPos = 0;
};

}  // namespace

/// \brief Detect if \p localRead is part of an anomalous read pair, and if so, convert the read pair
/// information into an SVObservation, and append it to \p candidates
///
/// \param[in] rstats Statistics computed from the fragment length distribution associated with \p localRead
///
/// \param[in] localRead The read which is the subject of testing/conversion to an SV candidate
///
/// \param[in] localAlign Pre-computed alignment data generated from \p localRead as a caching optimization
///
/// \param[in] remoteReadPtr Pointer to the bam record of \p localRead's mate. If nullptr, then properties of
/// the mate alignment are inferred from the local alignment record.
///
/// \param[in,out] candidates New SVObservation objects are appended to this vector. Contents of the vector
/// are preserved but not read.
static void getSVCandidatesFromPair(
    const ReadScannerOptions&                   opt,
    const ReadScannerDerivOptions&              dopt,
    const SVLocusScanner::CachedReadGroupStats& rstats,
    const bam_record&                           localRead,
    const SimpleAlignment&                      localAlign,
    const bam_record*                           remoteReadPtr,
    const bam_header_info&                      bamHeader,
    std::vector<SVObservation>&                 candidates)
{
  if (!localRead.is_paired()) return;

  // If this is one of the 'supplemental' sub-reads of a split read entry, don't use it here. This ensures
  // that each read pair is evaluated for anomolous pair evidence only once (when the non-supplemental read
  // is evaluated).
  if (localRead.isNonStrictSupplement()) return;

  // If either read is unmapped, an anomalous fragment size cannot be inferred.
  if (localRead.is_unmapped() || localRead.is_mate_unmapped()) return;

  // Special case typically used for RNA-Seq analysis:
  if (opt.isIgnoreAnomProperPair && localRead.is_proper_pair()) return;

  // All information about the remote alignment is consolidated to one object. If the remote alignment record
  // is available the remote alignment will be more accurate.
  const bool            isRemoteReadAvailable(nullptr != remoteReadPtr);
  const SimpleAlignment remoteAlign(
      isRemoteReadAvailable ? getAlignment(*remoteReadPtr) : getKnownOrFakedMateAlignment(localRead));

  AlignmentPairAnalyzer pairInspector(opt, dopt, rstats, bamHeader);
  pairInspector.reset(localAlign, remoteAlign, (!isRemoteReadAvailable), localRead.is_first());

  if (!pairInspector.isAnomalousReadPair()) return;

  if (pairInspector.isAlignedToChromEnds()) return;

  candidates.emplace_back();
  pairInspector.getSVObservation(candidates.back());

#ifdef DEBUG_SCANNER
  log_os << __FUNCTION__ << " evaluating pair sv for inclusion: " << candidates.back() << "\n";
#endif
}

#if 0
/// get SV candidates from shadow/singleton pairs
/// look for singletons, create candidateSV around conf. interval of shadow position
/// cache singletons? might be needed to remove poor quality shadows.
/// should be able to re-use code, follow soft-clipping example.
static
void
getSVCandidatesFromShadow(
    const ReadScannerOptions& opt,
    const SVLocusScanner::CachedReadGroupStats& rstats,
    const bam_record& localRead,
    const SimpleAlignment& localAlign,
    const bam_record* remoteReadPtr,
    TrackedCandidates& candidates)
{
    using namespace SVEvidenceType;
    static const index_t svSource(SHADOW);

    static const bool isComplex(true);
    pos_t singletonGenomePos(0);
    int targetId(0);
    if (nullptr == remoteReadPtr)
    {
        if (!localRead.is_unmapped()) return;
        // need to take care of this case
        // need to rely on cached mapq and qname
        return;
        if (!isGoodShadow(localRead,lastMapq,lastQname,opt.minSingletonMapqGraph))
        {
            return;
        }
        singletonGenomePos = localAlign.pos;
        targetId           = localRead.target_id();
    }
    else
    {
        // have both reads, straightforward from here
        const bam_record& remoteRead(*remoteReadPtr);
        const SimpleAlignment remoteAlign(remoteRead);

        if (localRead.is_mate_unmapped())
        {
            // remote read is shadow candidate
            if (!isGoodShadow(remoteRead,localRead.map_qual(),localRead.qname(),opt.minSingletonMapqGraph))
            {
                return;
            }
            singletonGenomePos = localAlign.pos;
            targetId = remoteRead.target_id();
        }
        else if (localRead.is_unmapped())
        {
            // local is shadow candidate
            if (!isGoodShadow(localRead,remoteRead.map_qual(),remoteRead.qname(),opt.minSingletonMapqGraph))
            {
                return;
            }
            singletonGenomePos = remoteAlign.pos;
            targetId = localRead.target_id();
        }
        else
        {
            // none unmapped, skip this one
            return;
        }
    }
    const pos_t properPairRangeOffset = static_cast<pos_t>(rstats.properPair.min + (rstats.properPair.max-rstats.properPair.min)/2);
    const pos_t shadowGenomePos = singletonGenomePos + properPairRangeOffset;
    candidates.push_back(GetSplitSVCandidate(opt,targetId,shadowGenomePos,shadowGenomePos, svSource, isComplex));
}
#endif

static void getSingleReadSVCandidates(
    const ReadScannerOptions&       opt,
    const ReadScannerDerivOptions&  dopt,
    const bam_record&               localRead,
    const SimpleAlignment&          localAlign,
    const bam_header_info&          bamHeader,
    const reference_contig_segment& refSeq,
    std::vector<SVObservation>&     candidates)
{
  using namespace illumina::common;

  const bool                                     isRead2(localRead.is_paired() && (localRead.read_no() == 2));
  const SourceOfSVEvidenceInDNAFragment::index_t fragSource(
      isRead2 ? SourceOfSVEvidenceInDNAFragment::READ2 : SourceOfSVEvidenceInDNAFragment::READ1);

  // - process any large indels in the localRead:
  getSVCandidatesFromReadIndels(opt, dopt, localAlign, bamHeader, fragSource, candidates);
#ifdef DEBUG_SCANNER
  log_os << __FUNCTION__ << ": post-indels candidate_size: " << candidates.size() << "\n";
#endif

  // a read can provide SA split evidence or semi-aligned/soft-clip, but not both.
  // this prevents split reads from triggering spurious local assembles. It is
  // possible for a read to genuinely contain evidence of both, but this should
  // be very rare.
  if (localRead.isSASplit()) {
    getSACandidatesFromRead(opt, dopt, localRead, localAlign, fragSource, bamHeader, candidates);
#ifdef DEBUG_SCANNER
    log_os << __FUNCTION__ << ": post-split read candidate_size: " << candidates.size() << "\n";
#endif
  } else {
    if (dopt.isSmallCandidates) {
      getSVCandidatesFromSemiAligned(
          opt, dopt, localRead, bamHeader, localAlign, fragSource, refSeq, candidates);
    }
#ifdef DEBUG_SCANNER
    log_os << __FUNCTION__ << ": post-semialigned candidate_size: " << candidates.size() << "\n";
#endif
  }
}

/// \brief Scan read record (and optionally its mate record) for SV evidence
//
/// Note that estimation is improved by the mate record (because we have the mate cigar string in this case)
///
static void getReadBreakendsImpl(
    const ReadScannerOptions&                   opt,
    const ReadScannerDerivOptions&              dopt,
    const SVLocusScanner::CachedReadGroupStats& rstats,
    const bam_record&                           localRead,
    const bam_record*                           remoteReadPtr,
    const bam_header_info&                      bamHeader,
    const reference_contig_segment&             localRefSeq,
    const reference_contig_segment*             remoteRefSeqPtr,
    std::vector<SVObservation>&                 candidates,
    known_pos_range2&                           localEvidenceRange)
{
  using namespace illumina::common;

#ifdef DEBUG_SCANNER
  log_os << __FUNCTION__ << ": Starting read: " << localRead.qname() << "\n";
#endif

  const chromMap_t& chromToIndex(bamHeader.chrom_to_index);

  candidates.clear();

  /// get some basic derived information from the bam_record:
  const SimpleAlignment localAlign(getAlignment(localRead));

  try {
    getSingleReadSVCandidates(opt, dopt, localRead, localAlign, bamHeader, localRefSeq, candidates);

    // run the same check on the read's mate if we have access to it
    if (nullptr != remoteReadPtr) {
      const bam_record&     remoteRead(*remoteReadPtr);
      const SimpleAlignment remoteAlign(getAlignment(remoteRead));

      if (nullptr == remoteRefSeqPtr) {
        static const char msg[] = "ERROR: remoteRefSeqPtr cannot be null";
        BOOST_THROW_EXCEPTION(GeneralException(msg));
      }
      getSingleReadSVCandidates(
          opt, dopt, remoteRead, remoteAlign, bamHeader, (*remoteRefSeqPtr), candidates);
    }

    // process shadows:
    //getSVCandidatesFromShadow(opt, rstats, localRead, localAlign,remoteReadPtr,candidates);

    // - process anomalous read pairs:
    getSVCandidatesFromPair(opt, dopt, rstats, localRead, localAlign, remoteReadPtr, bamHeader, candidates);
  } catch (...) {
    std::cerr << "ERROR: Exception caught while processing ";
    if (nullptr == remoteReadPtr) {
      std::cerr << "single read record:\n" << '\t' << localRead << "\n";
    } else {
      std::cerr << " read pair records:\n" << '\t' << localRead << "\n" << '\t' << (*remoteReadPtr) << "\n";
    }
    throw;
  }

#ifdef DEBUG_SCANNER
  log_os << __FUNCTION__ << ": post-pair candidate_size: " << candidates.size() << "\n";
#endif

  // update localEvidence range:
  // note this is only used if candidates were added, so there's no harm in setting it every time:
  const unsigned localRefLength(apath_ref_length(localAlign.path));
  const pos_t    startRefPos(localRead.pos() - 1);
  const pos_t    endRefPos(startRefPos + localRefLength);

  localEvidenceRange.set_range(startRefPos, endRefPos);

  const int maxTid(chromToIndex.size());

  /// final chance to QC candidate set:
  ///
  for (const SVCandidate& sv : candidates) {
    bool isInvalidTid(false);
    if ((sv.bp1.interval.tid < 0) || (sv.bp1.interval.tid >= maxTid)) {
      isInvalidTid = true;
    } else if (sv.bp2.state != SVBreakendState::UNKNOWN) {
      if ((sv.bp2.interval.tid < 0) || (sv.bp2.interval.tid >= maxTid)) {
        isInvalidTid = true;
      }
    }

    bool isInvalidPos(false);
    if (!isInvalidTid) {
      // note in the 'off-chromosome edge' test below we check for cases which are obviously way off
      // the edge, but allow for a bit of over-edge mistakes to occur for the circular chromosomes
      //
      static const int offEdgePad(500);
      const pos_t      tid1Length(bamHeader.chrom_data[sv.bp1.interval.tid].length);
      if ((sv.bp1.interval.range.end_pos() <= -offEdgePad) ||
          (sv.bp1.interval.range.begin_pos() >= (tid1Length + offEdgePad))) {
        isInvalidPos = true;
      } else if (sv.bp2.state != SVBreakendState::UNKNOWN) {
        const pos_t tid2Length(bamHeader.chrom_data[sv.bp2.interval.tid].length);
        if ((sv.bp2.interval.range.end_pos() <= -offEdgePad) ||
            (sv.bp2.interval.range.begin_pos() >= (tid2Length + offEdgePad))) {
          isInvalidPos = true;
        }
      }
    }

    if (isInvalidTid || isInvalidPos) {
      std::ostringstream oss;
      if (isInvalidTid) {
        oss << "SVbreakend has unknown or invalid chromosome id in candidate sv.\n";
      } else {
        oss << "Cannot interpret BAM record: candidate SV breakend from BAM record is off chromosome edge.\n";
      }

      oss << "\tlocal_bam_record: " << localRead << "\n"
          << "\tremote_bam record: ";
      if (nullptr == remoteReadPtr) {
        oss << "NONE";
      } else {
        oss << (*remoteReadPtr);
      }
      oss << "\n"
          << "\tSVCandidate: " << sv << "\n";
      BOOST_THROW_EXCEPTION(GeneralException(oss.str()));
    }
  }
}

/// \brief Create an SVLocus for each potential SV event supported by the BAM record
///
/// The loci count should almost always be one (or, depending on input filtration, zero).
/// multiple suggested loci from one read is more of a theoretical possibility than an
/// expectation.
///
static void getSVLociImpl(
    const ReadScannerOptions&                   opt,
    const ReadScannerDerivOptions&              dopt,
    const SVLocusScanner::CachedReadGroupStats& rstats,
    const bam_record&                           bamRead,
    const bam_header_info&                      bamHeader,
    const reference_contig_segment&             refSeq,
    std::vector<SVLocus>&                       loci,
    SampleEvidenceCounts&                       eCounts)
{
  using namespace illumina::common;

  loci.clear();
  std::vector<SVObservation> candidates;
  known_pos_range2           localEvidenceRange;

  getReadBreakendsImpl(
      opt, dopt, rstats, bamRead, nullptr, bamHeader, refSeq, nullptr, candidates, localEvidenceRange);

#ifdef DEBUG_SCANNER
  log_os << __FUNCTION__ << ": candidate_size: " << candidates.size() << "\n";
#endif

  // translate SVCandidate to a simpler form for use
  // in the SV locus graph:
  for (const SVCandidate& cand : candidates) {
    const bool isCandComplex(isComplexSV(cand));

    const SVBreakend& localBreakend(cand.bp1);
    const SVBreakend& remoteBreakend(cand.bp2);

    if ((0 == localBreakend.interval.range.size()) ||
        ((!isCandComplex) && (0 == remoteBreakend.interval.range.size()))) {
      std::ostringstream oss;
      oss << "Unexpected breakend pattern proposed from bam record.\n"
          << "\tlocal_breakend: " << localBreakend << "\n"
          << "\tremote_breakend: " << remoteBreakend << "\n"
          << "\tbam_record: " << bamRead << "\n";
      BOOST_THROW_EXCEPTION(GeneralException(oss.str()));
    }

    // update evidence stats:
    for (int i(0); i < SVEvidenceType::SIZE; ++i) {
      eCounts.eType[i] += localBreakend.lowresEvidence.getVal(i);
    }

    // determine the evidence weight of this candidate:
    unsigned localEvidenceWeight(0);
    unsigned remoteEvidenceWeight(0);

    if (localBreakend.getAnyNonPairCount() != 0) {
      localEvidenceWeight = SVObservationWeights::internalReadEvent;
      if (remoteBreakend.getAnyNonPairCount() != 0) {
        remoteEvidenceWeight = SVObservationWeights::internalReadEvent;
      }
    } else if (localBreakend.getLocalPairCount() != 0) {
      bool isClose(false);
      if (is_innie_pair(bamRead)) {
        isClose = (std::abs(bamRead.template_size()) < rstats.minDistantFragmentSize);
      }

      unsigned thisWeight(SVObservationWeights::readPair);
      if (isClose) {
        thisWeight = SVObservationWeights::closeReadPair;
        eCounts.closeCount += 1;
      }

      localEvidenceWeight = thisWeight;
      if (remoteBreakend.getLocalPairCount() != 0) {
        remoteEvidenceWeight = thisWeight;
      }
    }

    // finally, create the graph locus:
    SVLocus locus;
    // set local breakend estimate:
    const NodeIndexType localBreakendNode(locus.addNode(localBreakend.interval));
    locus.setNodeEvidence(localBreakendNode, localEvidenceRange);

    if (isCandComplex) {
      locus.linkNodes(localBreakendNode, localBreakendNode, localEvidenceWeight);
    } else {
      // set remote breakend estimate:
      const NodeIndexType remoteBreakendNode(locus.addNode(remoteBreakend.interval));
      locus.linkNodes(localBreakendNode, remoteBreakendNode, localEvidenceWeight, remoteEvidenceWeight);

      locus.mergeSelfOverlap();
    }

#ifdef DEBUG_SCANNER
    log_os << __FUNCTION__ << ": adding Locus: " << locus << "\n";
#endif
    loci.push_back(locus);
  }
}

/// \brief Given a size distribution and probability \p prob, set \p range to span quantile(prob) to
/// quantile(1-prob)
///
/// \param[in] prob Probability used to determine quantile range drawn from the \p sizeDistribution
/// \param[out] range Computed quantile range over size
static void getSizeQuantileRange(
    const SizeDistribution& sizeDistribution, const float prob, SVLocusScanner::Range& range)
{
  range.min = sizeDistribution.quantile(prob);
  range.max = sizeDistribution.quantile((1 - prob));
  if (range.min < 0.) range.min = 0;
  assert(range.max > 0.);
}

SVLocusScanner::SVLocusScanner(
    const ReadScannerOptions& opt,
    const std::string&        statsFilename,
    const std::vector<std::string>& /*alignmentFilename*/,
    const bool isTranscriptStrandKnown)
  : _opt(opt), _dopt(opt, isTranscriptStrandKnown)
{
  using namespace illumina::common;

  // pull in insert stats:
  _rss.load(statsFilename.c_str());

  // precompute frequently used insert stats for each read group:
  const unsigned readGroupCount(_rss.size());
  for (unsigned readGroupIndex(0); readGroupIndex < readGroupCount; readGroupIndex++) {
    /// \TODO add check that the filenames in the stats file are a complete match to alignmentFilename

    const SizeDistribution& readGroupFragSizeDistro(getFragSizeDistro(readGroupIndex));

    _stats.resize(_stats.size() + 1);
    CachedReadGroupStats& rgStats(_stats.back());
    getSizeQuantileRange(readGroupFragSizeDistro, _opt.breakendEdgeQuantileProb, rgStats.breakendRegion);
    getSizeQuantileRange(
        readGroupFragSizeDistro,
        _opt.largeScaleEventBreakendEdgeQuantileProb,
        rgStats.largeScaleEventBreakendRegion);
    getSizeQuantileRange(readGroupFragSizeDistro, _opt.properPairQuantileProb, rgStats.properPair);
    getSizeQuantileRange(readGroupFragSizeDistro, _opt.evidenceTrimQuantileProb, rgStats.evidencePair);
    getSizeQuantileRange(readGroupFragSizeDistro, 0.05f, rgStats.fifthPerc);

    if ((readGroupIndex == 0) || (rgStats.fifthPerc.min < _fifthPerc.min)) {
      _fifthPerc.min = rgStats.fifthPerc.min;
    }
    if ((readGroupIndex == 0) || (rgStats.fifthPerc.max > _fifthPerc.max)) {
      _fifthPerc.max = rgStats.fifthPerc.max;
    }

    rgStats.shadowSearchDistance =
        readGroupFragSizeDistro.quantile(1 - (_opt.shadowSearchDistanceQuantileProb)) *
        _opt.shadowSearchDistanceFactor;

    assert(rgStats.shadowSearchDistance > 0);

    rgStats.minDistantFragmentSize =
        static_cast<int>(rgStats.properPair.max * FragmentSizeType::closePairFactor);

    assert(rgStats.minDistantFragmentSize > rgStats.properPair.max);

    const int largeEventRegionMin(rgStats.properPair.max * FragmentSizeType::minLargeEventRegionFactor);
    const int largeEventRegionMax(rgStats.properPair.max * FragmentSizeType::maxLargeEventRegionFactor);

    rgStats.largeEventRegionScaler.init(largeEventRegionMin, largeEventRegionMax);
  }
}

void SVLocusScanner::checkReadSize(const stream_state_reporter& alignmentStream, const bam_record& bamRead)
{
  ALIGNPATH::path_t path;
  bam_cigar_to_apath(bamRead.raw_cigar(), bamRead.n_cigar(), path);
  const unsigned alignedSize(ALIGNPATH::apath_read_length(path));
  const unsigned seqSize(bamRead.read_size());

  if (seqSize == 0) {
    std::ostringstream oss;
    oss << "Input alignment record contains unknown read sequence (SEQ='*'), "
        << "which cannot be used for variant calling:\n";
    alignmentStream.report_state(oss);
    oss << "\n";
    BOOST_THROW_EXCEPTION(illumina::common::GeneralException(oss.str()));
  }

  if (seqSize != alignedSize) {
    std::ostringstream oss;
    oss << "Read length implied by mapped alignment (" << alignedSize << ") does not match sequence length ("
        << seqSize << "):\n";
    alignmentStream.report_state(oss);
    BOOST_THROW_EXCEPTION(illumina::common::GeneralException(oss.str()));
  }
}

bool SVLocusScanner::isAnomalousReadPair(
    const bam_record& bamRead, const unsigned defaultReadGroupIndex) const
{
  if (!is_innie_pair(bamRead)) return true;

  const Range&  ppr(_stats[defaultReadGroupIndex].properPair);
  const int32_t fragmentSize(std::abs(bamRead.template_size()));

  // Unknown fragment size
  if (fragmentSize == 0) return true;

  // we're seeing way to much large fragment garbage in cancers to use
  // vanilla proper pair criteria, push the max fragment size out a bit for now:
  static const float maxAnomFactor(1.5);
  if (fragmentSize > static_cast<int32_t>(maxAnomFactor * ppr.max)) return true;
  if (fragmentSize < ppr.min) return true;

  return false;
}

FragmentSizeType::index_t SVLocusScanner::_getFragmentSizeType(
    const bam_record& bamRead, const unsigned defaultReadGroupIndex) const
{
  using namespace FragmentSizeType;
  if (bamRead.target_id() != bamRead.mate_target_id()) return DISTANT;
  const int32_t fragmentSize(std::abs(bamRead.template_size()));
  if (fragmentSize == 0) return UNKNOWN;
  return classifySize(_stats[defaultReadGroupIndex], fragmentSize);
}

bool SVLocusScanner::_isLargeFragment(const bam_record& bamRead, const unsigned defaultReadGroupIndex) const
{
  return FragmentSizeType::isLarge(_getFragmentSizeType(bamRead, defaultReadGroupIndex));
}

bool SVLocusScanner::isNonCompressedAnomalousReadPair(
    const bam_record& bamRead, const unsigned defaultReadGroupIndex) const
{
  if (!is_mapped_pair(bamRead)) return false;
  const bool isAnomalous(isAnomalousReadPair(bamRead, defaultReadGroupIndex));
  const bool isInnie(is_innie_pair(bamRead));
  const bool isLarge(_isLargeFragment(bamRead, defaultReadGroupIndex));

  // exclude innie read pairs which are anomalously short:
  return (isAnomalous && ((!isInnie) || isLarge));
}

bool SVLocusScanner::isLocalIndelEvidence(const SimpleAlignment& bamAlign) const
{
  using namespace ALIGNPATH;
  for (const path_segment& ps : bamAlign.path) {
    if (ps.type == INSERT || ps.type == DELETE) {
      if (ps.length >= _opt.minCandidateVariantSize) return true;
    }
  }
  return false;
}

bool SVLocusScanner::isSemiAlignedEvidence(
    const bam_record& bamRead, const SimpleAlignment& bamAlign, const reference_contig_segment& refSeq) const
{
  unsigned leadingMismatchLen(0), trailingMismatchLen(0);
  getSVBreakendCandidateSemiAlignedSimple(
      bamRead, bamAlign, refSeq, _opt.useOverlapPairEvidence, leadingMismatchLen, trailingMismatchLen);
  return (
      (leadingMismatchLen >= _opt.minSemiAlignedMismatchLen) ||
      (trailingMismatchLen >= _opt.minSemiAlignedMismatchLen));
}

bool SVLocusScanner::isLocalAssemblyEvidence(
    const bam_record& bamRead, const reference_contig_segment& refSeq) const
{
  const SimpleAlignment bamAlign(getAlignment(bamRead));
  if (isLocalIndelEvidence(bamAlign)) return true;
  if (isSemiAlignedEvidence(bamRead, bamAlign, refSeq)) return true;
  /// TODO Add shadow evidence -- complexity here is keeping locus merging under control due to the large
  /// breakend location variance suggested by shadows

  return false;
}

bool SVLocusScanner::isSVEvidence(
    const bam_record&               bamRead,
    const unsigned                  defaultReadGroupIndex,
    const reference_contig_segment& refSeq,
    SVLocusEvidenceCount*           incountsPtr) const
{
  // exclude innie read pairs which are anomalously short:
  const bool      isAnom(isNonCompressedAnomalousReadPair(bamRead, defaultReadGroupIndex));
  const bool      isSplit(bamRead.isSASplit());
  SimpleAlignment _bamAlign;
  getAlignment(bamRead, _bamAlign);
  const bool isIndel(isLocalIndelEvidence(_bamAlign));
  const bool isAssm(
      (_dopt.isSmallCandidates) && ((!isSplit) && isSemiAlignedEvidence(bamRead, _bamAlign, refSeq)));

  // Mark supplemental segments of split reads, and exclude these from the "normal" read counts
  const bool isSupplementary(bamRead.is_supplementary() || bamRead.is_secondary());

  const bool isEvidence(isAnom || isSplit || isIndel || isAssm);

  if (nullptr != incountsPtr) {
    SVLocusEvidenceCount& incounts(*incountsPtr);

    if (isSupplementary) {
      assert(isSplit);
      incounts.splitSupplementarySegment++;
    } else {
      incounts.total++;
      if (isAnom) incounts.anom++;
      if (isSplit) incounts.split++;
      if (isAnom && isSplit) incounts.anomAndSplit++;
      if (isIndel) incounts.indel++;
      if (isAssm) incounts.assm++;

      if (!isEvidence) incounts.ignored++;

      if (isAnom) {
        if (isMateInsertionEvidenceCandidate(bamRead, getMinMapQ())) {
          // these counts are used to generate background noise rates in later candidate generation stages:
          incounts.remoteRecoveryCandidates++;
        }
      }
    }
  }

  return isEvidence;
}

void SVLocusScanner::getSVLoci(
    const bam_record&               bamRead,
    const unsigned                  defaultReadGroupIndex,
    const bam_header_info&          bamHeader,
    const reference_contig_segment& refSeq,
    std::vector<SVLocus>&           loci,
    SampleEvidenceCounts&           eCounts) const
{
  loci.clear();

  const CachedReadGroupStats& rstats(_stats[defaultReadGroupIndex]);
  getSVLociImpl(_opt, _dopt, rstats, bamRead, bamHeader, refSeq, loci, eCounts);
}

void SVLocusScanner::getBreakendPair(
    const bam_record&               localRead,
    const bam_record*               remoteReadPtr,
    const unsigned                  defaultReadGroupIndex,
    const bam_header_info&          bamHeader,
    const reference_contig_segment& localRefSeq,
    const reference_contig_segment* remoteRefSeqPtr,
    std::vector<SVObservation>&     candidates) const
{
  const CachedReadGroupStats& rstats(_stats[defaultReadGroupIndex]);

  // throw evidence range away in this case
  known_pos_range2 evidenceRange;
  getReadBreakendsImpl(
      _opt,
      _dopt,
      rstats,
      localRead,
      remoteReadPtr,
      bamHeader,
      localRefSeq,
      remoteRefSeqPtr,
      candidates,
      evidenceRange);
}
