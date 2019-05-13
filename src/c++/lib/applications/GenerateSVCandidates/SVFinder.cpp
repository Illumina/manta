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
/// \author Naoki Nariai
///

#include "SVFinder.hpp"

#include <iostream>

#include "blt_util/binomial_test.hpp"
#include "blt_util/log.hpp"
#include "common/Exceptions.hpp"
#include "manta/BamStreamerUtils.hpp"
#include "manta/ReadFilter.hpp"
#include "manta/SVCandidateUtil.hpp"
#include "manta/SVReferenceUtil.hpp"
#include "svgraph/EdgeInfoUtil.hpp"

//#define DEBUG_SVDATA

#ifdef DEBUG_SVDATA
#include "blt_util/log.hpp"
#endif

static double getSpanningNoiseRate(const AllSampleReadCounts& counts, const unsigned sampleIndex)
{
  static const double pseudoTotal(1000.);
  static const double pseudoSpan(10.);

  const SampleReadInputCounts& input(counts.getSampleCounts(sampleIndex).input);

  const double anomOrSplitCount(
      input.evidenceCount.anom + input.evidenceCount.split - input.evidenceCount.anomAndSplit);
  assert(anomOrSplitCount >= 0);
  assert(anomOrSplitCount <= input.evidenceCount.total);

  return (anomOrSplitCount + pseudoSpan) / (input.evidenceCount.total + pseudoTotal);
}

static double getAssemblyNoiseRate(const AllSampleReadCounts& counts, const unsigned sampleIndex)
{
  static const double pseudoTotal(1000.);
  static const double pseudoAssm(10.);

  const SampleReadInputCounts& input(counts.getSampleCounts(sampleIndex).input);
  assert(input.evidenceCount.assm <= input.evidenceCount.total);

  return (input.evidenceCount.assm + pseudoAssm) / (input.evidenceCount.total + pseudoTotal);
}

SVFinder::SVFinder(
    const GSCOptions&                   opt,
    const SVLocusScanner&               readScanner,
    const bam_header_info&              bamHeader,
    const AllSampleReadCounts&          readCounts,
    std::shared_ptr<EdgeRuntimeTracker> edgeTrackerPtr,
    GSCEdgeStatsManager&                edgeStatMan)
  : _scanOpt(opt.scanOpt),
    _isAlignmentTumor(opt.alignFileOpt.isAlignmentTumor),
    _readScanner(readScanner),
    _referenceFilename(opt.referenceFilename),
    _skipEvidenceSignalFilter(opt.skipEvidenceSignalFilter),
    _isRNA(opt.isRNA),
    _isVerbose(opt.isVerbose),
    _isSomatic(false),
    _edgeTrackerPtr(edgeTrackerPtr),
    _edgeStatMan(edgeStatMan)
{
  _dFilterPtr.reset(new ChromDepthFilterUtil(opt.chromDepthFilename, _scanOpt.maxDepthFactor, bamHeader));

  // setup regionless bam_streams:
  openBamStreams(opt.referenceFilename, opt.alignFileOpt.alignmentFilenames, _bamStreams);

  const unsigned bamCount(_bamStreams.size());
  {
    // assert expected bam order of all normal samples followed by all tumor samples,
    // also, determine if this is a somatic run:
    bool isFirstTumor(false);
    for (unsigned bamIndex(0); bamIndex < bamCount; ++bamIndex) {
      const bool isTumor(_isAlignmentTumor[bamIndex]);

      if (isTumor) {
        isFirstTumor = true;
        _isSomatic   = true;
      }
      assert((!isFirstTumor) || isTumor);
    }
  }

  for (unsigned bamIndex(0); bamIndex < bamCount; ++bamIndex) {
    _spanningNoiseRate.push_back(getSpanningNoiseRate(readCounts, bamIndex));
    _assemblyNoiseRate.push_back(getAssemblyNoiseRate(readCounts, bamIndex));
  }
}

// making the dtor explicit and in the cpp allows unique_ptr to work reliably:
SVFinder::~SVFinder() {}

/// test if read supports an SV on this edge, if so, add to SVData
static void addSVNodeRead(
    const bam_header_info&                     bamHeader,
    const SVLocusScanner&                      scanner,
    const SVLocusNode&                         localNode,
    const SVLocusNode&                         remoteNode,
    const bam_record&                          bamRead,
    const unsigned                             bamIndex,
    const bool                                 isExpectRepeat,
    const reference_contig_segment&            refSeq,
    const bool                                 isLocalNodeGraphEdgeNode1,
    const bool                                 isGatherSubmapped,
    SVCandidateSetSequenceFragmentSampleGroup& svDataGroup,
    SampleEvidenceCounts&                      eCounts)
{
  using namespace illumina::common;

  if (bamRead.map_qual() < scanner.getMinTier2MapQ()) return;

  const bool isSubMapped(bamRead.map_qual() < scanner.getMinMapQ());
  if ((!isGatherSubmapped) && isSubMapped) return;

  svDataGroup.increment(isSubMapped);

  if (!scanner.isSVEvidence(bamRead, bamIndex, refSeq)) return;

  // finally, check to see if the svDataGroup is full... for now, we allow a very large
  // number of reads to be stored in the hope that we never reach this limit, but just in
  // case we don't want to exhaust memory in centromere pileups, etc...
  //
  // Once svDataGroup is full, we keep scanning but only to find pairs for reads that
  // have already been entered.
  //
  static const unsigned maxDataSize(4000);
  if ((!svDataGroup.isFull()) && (svDataGroup.size() >= maxDataSize)) {
    svDataGroup.setFull();
  }

  //
  // run an initial screen to make sure at least one candidate from this read matches the regions for this
  // edge:
  //
  typedef std::vector<SVLocus> loci_t;
  loci_t                       loci;
  scanner.getSVLoci(bamRead, bamIndex, bamHeader, refSeq, loci, eCounts);

  for (const SVLocus& locus : loci) {
    const unsigned locusSize(locus.size());
    assert((locusSize >= 1) && (locusSize <= 2));

    unsigned readLocalIndex(0);
    if (locusSize == 2) {
      unsigned readRemoteIndex(1);
      if (!locus.getNode(readLocalIndex).isOutCount()) {
        std::swap(readLocalIndex, readRemoteIndex);
      }

      if (!locus.getNode(readLocalIndex).isOutCount()) {
        std::ostringstream oss;
        oss << "Unexpected svlocus counts from bam record: " << bamRead << "\n"
            << "\tlocus: " << locus << "\n";
        BOOST_THROW_EXCEPTION(GeneralException(oss.str()));
      }
      if (!locus.getNode(readRemoteIndex).getInterval().isIntersect(remoteNode.getInterval()))
        continue;  // todo should this intersect be checked in swapped orientation?
    } else {
      if (!locus.getNode(readLocalIndex).getInterval().isIntersect(remoteNode.getInterval())) continue;
    }

    if (!locus.getNode(readLocalIndex).getInterval().isIntersect(localNode.getInterval()))
      continue;  // todo should this intersect be checked in swapped orientation?

    svDataGroup.add(bamHeader, bamRead, isExpectRepeat, isLocalNodeGraphEdgeNode1, isSubMapped);

    // once any loci has achieved the local/remote overlap criteria, there's no reason to keep scanning loci
    // of the same bam record:
    break;
  }
}

/// \brief Load a segment of the reference encompassing a node from the SV locus graph
static void getNodeRefSeq(
    const bam_header_info&    bamHeader,
    const SVLocus&            locus,
    const NodeIndexType       localNodeIndex,
    const std::string&        referenceFilename,
    GenomeInterval&           searchInterval,
    reference_contig_segment& refSeq)
{
  // get full search interval:
  const SVLocusNode& localNode(locus.getNode(localNodeIndex));
  searchInterval = (localNode.getInterval());
  searchInterval.range.merge_range(localNode.getEvidenceRange());

  // grab the reference for segment we're estimating plus a buffer around the segment edges:
  static const unsigned refEdgeBufferSize(100);
  getIntervalReferenceSegment(referenceFilename, bamHeader, refEdgeBufferSize, searchInterval, refSeq);
}

/// approximate depth tracking -- don't bother reading the cigar string, just assume a perfect match of
/// size read_size
static void addReadToDepthEst(const bam_record& bamRead, const pos_t beginPos, std::vector<unsigned>& depth)
{
  const pos_t endPos(beginPos + depth.size());
  const pos_t refStart(bamRead.pos() - 1);

  const pos_t readSize(bamRead.read_size());
  for (pos_t readIndex(std::max(0, (beginPos - refStart))); readIndex < readSize; ++readIndex) {
    const pos_t refPos(refStart + readIndex);
    if (refPos >= endPos) return;
    const pos_t depthIndex(refPos - beginPos);
    assert(depthIndex >= 0);

    depth[depthIndex]++;
  }
}

void SVFinder::addSVNodeData(
    const bam_header_info&          bamHeader,
    const SVLocus&                  locus,
    const NodeIndexType             localNodeIndex,
    const NodeIndexType             remoteNodeIndex,
    const GenomeInterval&           searchInterval,
    const reference_contig_segment& localNodeRefSeq,
    const bool                      isLocalNodeGraphEdgeNode1,
    SVCandidateSetData&             svData)
{
  // get full search interval:
  const SVLocusNode& localNode(locus.getNode(localNodeIndex));
  const SVLocusNode& remoteNode(locus.getNode(remoteNodeIndex));

  bool isExpectRepeat(svData.setNewSearchInterval(searchInterval));

  // This is a temporary measure to make the QNAME collision detection much less strict.
  //
  // Problems have come up where very large deletions are present in a read, and it is therefore
  // detected as a repeat in two different regions, even though they are separated by a considerable
  // distance. The solution below is to temporarily turn off QNAME collision detection whenever two
  // regions are on the same chrom (ie. almost always)
  //
  // TODO restore more precise collision detection
  if (!isExpectRepeat) isExpectRepeat = (localNode.getInterval().tid == remoteNode.getInterval().tid);

#ifdef DEBUG_SVDATA
  log_os << __FUNCTION__ << ": bp_interval: " << localNode.getInterval()
         << " evidenceInterval: " << localNode.getEvidenceRange() << " searchInterval: " << searchInterval
         << " isExpectRepeat: " << isExpectRepeat << "\n";
#endif

  const bool isMaxDepth(dFilter().isMaxDepthFilter());
  float      maxDepth(0);
  if (isMaxDepth) {
    maxDepth = dFilter().maxDepth(searchInterval.tid);
  }
  const pos_t           searchBeginPos(searchInterval.range.begin_pos());
  const pos_t           searchEndPos(searchInterval.range.end_pos());
  std::vector<unsigned> normalDepthBuffer(searchInterval.range.size(), 0);

  // iterate through reads, test reads for association and add to svData:
  unsigned bamIndex(0);
  for (streamPtr& bamPtr : _bamStreams) {
    const bool isTumor(_isAlignmentTumor[bamIndex]);

    const bool isGatherSubmapped(_isSomatic && (!isTumor));

    SVCandidateSetSequenceFragmentSampleGroup& svDataGroup(svData.getDataGroup(bamIndex));
    bam_streamer&                              readStream(*bamPtr);

    // record the associated stream name to improve error messages
    svDataGroup.dataSourceName = readStream.name();

    // set bam stream to new search interval:
    readStream.resetRegion(
        searchInterval.tid, searchInterval.range.begin_pos(), searchInterval.range.end_pos());

#ifdef DEBUG_SVDATA
    log_os << __FUNCTION__ << ": scanning bamIndex: " << bamIndex << "\n";
#endif
    while (readStream.next()) {
      const bam_record& bamRead(*(readStream.get_record_ptr()));

      const pos_t refPos(bamRead.pos() - 1);
      if (refPos >= searchEndPos) break;

      if (isReadUnmappedOrFilteredCore(bamRead)) continue;

      if (isMaxDepth) {
        if (!isTumor) {
          addReadToDepthEst(bamRead, searchBeginPos, normalDepthBuffer);
        }

        assert(refPos < searchEndPos);
        const pos_t depthOffset(refPos - searchBeginPos);
        if ((depthOffset >= 0) && (normalDepthBuffer[depthOffset] > maxDepth)) continue;
      }

      // test if read supports an SV on this edge, if so, add to SVData
      addSVNodeRead(
          bamHeader,
          _readScanner,
          localNode,
          remoteNode,
          bamRead,
          bamIndex,
          isExpectRepeat,
          localNodeRefSeq,
          isLocalNodeGraphEdgeNode1,
          isGatherSubmapped,
          svDataGroup,
          _eCounts);
    }
    ++bamIndex;
  }
}

// sanity check the final result
void SVFinder::checkResult(const SVCandidateSetData& svData, const std::vector<SVCandidate>& svs) const
{
  using namespace illumina::common;

  const unsigned svCount(svs.size());
  if (0 == svCount) return;

  // check that the counts totaled up from the data match those in the sv candidates
  std::map<unsigned, unsigned> readCounts;
  std::map<unsigned, unsigned> pairCounts;

  for (unsigned i(0); i < svCount; ++i) {
    readCounts[i] = 0;
    pairCounts[i] = 0;
  }

  const unsigned bamCount(_bamStreams.size());
  for (unsigned bamIndex(0); bamIndex < bamCount; ++bamIndex) {
    const SVCandidateSetSequenceFragmentSampleGroup& svDataGroup(svData.getDataGroup(bamIndex));
    for (const SVCandidateSetSequenceFragment& fragment : svDataGroup) {
      for (const SVSequenceFragmentAssociation& sva : fragment.svLink) {
        if (sva.index >= svCount) {
          std::ostringstream oss;
          oss << "Searching for SVIndex: " << sva.index << " with svSize: " << svCount << "\n";
          BOOST_THROW_EXCEPTION(GeneralException(oss.str()));
        }

        if (SVEvidenceType::isPairType(sva.evtype)) {
          if (fragment.read1.isSet()) readCounts[sva.index]++;
          if (fragment.read2.isSet()) readCounts[sva.index]++;
          if (fragment.read1.isSet() && fragment.read2.isSet()) pairCounts[sva.index] += 2;
        }
      }
    }
  }

  for (unsigned svIndex(0); svIndex < svCount; ++svIndex) {
    const unsigned svObsReadCount(
        svs[svIndex].bp1.getLocalPairCount() + svs[svIndex].bp2.getLocalPairCount());
    const unsigned svObsPairCount(svs[svIndex].bp1.getPairCount() + svs[svIndex].bp2.getPairCount());
    assert(svs[svIndex].bp1.getPairCount() == svs[svIndex].bp2.getPairCount());

    const unsigned dataObsReadCount(readCounts[svIndex]);
    const unsigned dataObsPairCount(pairCounts[svIndex]);

    bool isCountException(false);
    if (svObsReadCount != dataObsReadCount)
      isCountException = true;
    else if (svObsPairCount != dataObsPairCount)
      isCountException = true;

    if (isCountException) {
      std::ostringstream oss;
      oss << "Unexpected difference in sv and data read counts.\n"
          << "\tSVreadCount: " << svObsReadCount << " DataReadCount: " << dataObsReadCount << "\n"
          << "\tSVpaircount: " << svObsPairCount << " DataPaircount: " << dataObsPairCount << "\n"
          << "\tsvIndex: " << svIndex << " SV: " << svs[svIndex];
      BOOST_THROW_EXCEPTION(GeneralException(oss.str()));
    }
  }
}

typedef std::map<unsigned, unsigned> movemap_t;

/// local convenience struct, if only I had closures instead... :<
struct svCandDeleter {
  svCandDeleter(std::vector<FatSVCandidate>& svs, movemap_t& moveSVIndex)
    : _shift(0), _isLastIndex(false), _lastIndex(0), _svs(svs), _moveSVIndex(moveSVIndex)
  {
  }

  void deleteIndex(const unsigned index)
  {
    assert(index <= _svs.size());

    if (_isLastIndex) {
      for (unsigned i(_lastIndex + 1); i < index; ++i) {
        assert(_shift > 0);
        assert(i >= _shift);

        _svs[(i - _shift)] = _svs[i];
        // moveSVIndex has already been set for deleted indices, this sets
        // the move for non-deleted positions:
        _moveSVIndex[i] = (i - _shift);
      }
    }
    _lastIndex   = index;
    _isLastIndex = true;
    _shift++;
  }

private:
  unsigned                     _shift;
  bool                         _isLastIndex;
  unsigned                     _lastIndex;
  std::vector<FatSVCandidate>& _svs;
  movemap_t&                   _moveSVIndex;
};

/// check whether any svs have grown to intersect each other
///
/// this is also part of the temp hygen hack, so just make this minimally work:
///
static void consolidateOverlap(
    const unsigned bamCount, SVCandidateSetData& svData, std::vector<FatSVCandidate>& svs)
{
  movemap_t          moveSVIndex;
  std::set<unsigned> deletedSVIndex;

  std::vector<unsigned> innerIndexShift;

  const unsigned svCount(svs.size());
  for (unsigned outerIndex(1); outerIndex < svCount; ++outerIndex) {
    const unsigned prevInnerIndexShift((outerIndex <= 1) ? 0 : innerIndexShift[outerIndex - 2]);
    innerIndexShift.push_back(prevInnerIndexShift + deletedSVIndex.count(outerIndex - 1));
    for (unsigned innerIndex(0); innerIndex < outerIndex; ++innerIndex) {
      if (deletedSVIndex.count(innerIndex)) continue;

      if (svs[innerIndex].isIntersect(svs[outerIndex])) {
#ifdef DEBUG_SVDATA
        log_os << __FUNCTION__ << ": Merging outer:inner: " << outerIndex << " " << innerIndex << "\n";
#endif
        svs[innerIndex].merge(svs[outerIndex]);
        assert(innerIndexShift.size() > innerIndex);
        assert(innerIndexShift[innerIndex] <= innerIndex);
        moveSVIndex[outerIndex] = (innerIndex - innerIndexShift[innerIndex]);
        deletedSVIndex.insert(outerIndex);
        break;
      }
    }
  }

  if (!deletedSVIndex.empty()) {
#ifdef DEBUG_SVDATA
    for (const unsigned index : deletedSVIndex) {
      log_os << __FUNCTION__ << ": deleted index: " << index << "\n";
    }
#endif

    {
      svCandDeleter svDeleter(svs, moveSVIndex);

      for (const unsigned index : deletedSVIndex) {
        svDeleter.deleteIndex(index);
      }
      svDeleter.deleteIndex(svCount);
    }

    svs.resize(svs.size() - deletedSVIndex.size());

    // fix indices:
    for (unsigned i(0); i < svs.size(); ++i) {
      svs[i].candidateIndex = i;
    }
  }

  if (!moveSVIndex.empty()) {
#ifdef DEBUG_SVDATA
    for (const movemap_t::value_type& val : moveSVIndex) {
      log_os << __FUNCTION__ << ": Movemap from: " << val.first << " to: " << val.second << "\n";
    }
#endif

    for (unsigned bamIndex(0); bamIndex < bamCount; ++bamIndex) {
      SVCandidateSetSequenceFragmentSampleGroup& svDataGroup(svData.getDataGroup(bamIndex));
      for (SVCandidateSetSequenceFragment& fragment : svDataGroup) {
        for (SVSequenceFragmentAssociation& sva : fragment.svLink) {
          if (moveSVIndex.count(sva.index)) {
            sva.index = moveSVIndex[sva.index];
          }
        }
      }
    }
  }
}

/// Store additional signal rate information to help decide if the candidate evidence
/// is significant relative to background noise in the sample.
///
/// TODO -- is the following comment out of date?: only handles complex cases for now (assumes sv is complex)
static void updateEvidenceIndex(
    const SVCandidateSetSequenceFragment& fragment,
    const SVObservation&                  obs,
    FatSVCandidate&                       sv,
    const unsigned                        bamIndex)
{
  if (obs.isSingleReadSource()) {
    const SVCandidateSetRead& candRead(obs.isRead1Source() ? fragment.read1 : fragment.read2);
    if (obs.svEvidenceType != SVEvidenceType::SPLIT_ALIGN) {
      if ((not candRead.bamrec.empty()) && (not candRead.isSubMapped)) {
        sv.bp1EvidenceIndex[obs.svEvidenceType][bamIndex].push_back(candRead.readIndex);
#ifdef DEBUG_SVDATA
        log_os << __FUNCTION__ << ": non_split: Added readIndex " << candRead.readIndex
               << " to svEvidenceType " << obs.svEvidenceType << "\n";
#endif
      }
    } else {
      const bool  is1to1(sv.isIntersect1to1(obs));
      auto&       readBp(is1to1 ? sv.bp1EvidenceIndex : sv.bp2EvidenceIndex);
      auto&       readSuppBp(is1to1 ? sv.bp2EvidenceIndex : sv.bp1EvidenceIndex);
      const auto& read(obs.isRead1Source() ? fragment.read1 : fragment.read2);
      const auto& readSupp(obs.isRead1Source() ? fragment.read1Supplemental : fragment.read2Supplemental);

      if ((not read.bamrec.empty()) && (not read.isSubMapped)) {
        readBp[obs.svEvidenceType][bamIndex].push_back(read.readIndex);
#ifdef DEBUG_SVDATA
        log_os << __FUNCTION__ << ": split_mapped: Added readIndex " << candRead.readIndex
               << " to svEvidenceType " << obs.svEvidenceType << "\n";
#endif
      }
      if (readSupp.size() == 1) {
        if ((not readSupp.front().bamrec.empty()) && (not readSupp.front().isSubMapped)) {
          readSuppBp[obs.svEvidenceType][bamIndex].push_back(readSupp.front().readIndex);
#ifdef DEBUG_SVDATA
          log_os << __FUNCTION__ << ": split_supp_mapped: Added readIndex " << readSupp.front().readIndex
                 << " to svEvidenceType " << obs.svEvidenceType << "\n";
#endif
        }
      }
    }
  } else {
    // account for bp1 and bp2 mapping to read1 and read2:
    const bool                is1to1(sv.isIntersect1to1(obs));
    const SVCandidateSetRead& bp1Read(is1to1 ? fragment.read1 : fragment.read2);
    const SVCandidateSetRead& bp2Read(is1to1 ? fragment.read2 : fragment.read1);
    if ((not bp1Read.bamrec.empty()) && (not bp1Read.isSubMapped)) {
      sv.bp1EvidenceIndex[obs.svEvidenceType][bamIndex].push_back(bp1Read.readIndex);
#ifdef DEBUG_SVDATA
      log_os << __FUNCTION__ << ": multi_read_source: Added readIndex " << bp1Read.readIndex
             << " to svEvidenceType " << obs.svEvidenceType << "\n";
#endif
    }
    if ((not bp1Read.bamrec.empty()) && (not bp2Read.isSubMapped)) {
      sv.bp2EvidenceIndex[obs.svEvidenceType][bamIndex].push_back(bp2Read.readIndex);
#ifdef DEBUG_SVDATA
      log_os << __FUNCTION__ << ": multi_read_source: Added readIndex " << bp2Read.readIndex
             << " to svEvidenceType " << obs.svEvidenceType << "\n";
#endif
    }
  }
}

void SVFinder::assignFragmentObservationsToSVCandidates(
    const SVLocusNode&                node1,
    const SVLocusNode&                node2,
    const std::vector<SVObservation>& readCandidates,
    const bool                        isExpandSVCandidateSet,
    SVCandidateSetSequenceFragment&   fragment,
    std::vector<FatSVCandidate>&      svs,
    const unsigned                    bamIndex)
{
  const unsigned bamCount(_bamStreams.size());

  // we anticipate so few svs from the POC method, that there's no indexing on them
  for (const SVObservation& readCand : readCandidates) {
#ifdef DEBUG_SVDATA
    log_os << __FUNCTION__ << ": Starting assignment for read cand: " << readCand << "\n";
#endif
    if (_isRNA) {
      int minLength = isCis(readCand) ? _scanOpt.minRNACisLength : _scanOpt.minRNALength;
      if (isSVBelowMinSize(readCand, minLength)) {
#ifdef DEBUG_SVDATA
        log_os << __FUNCTION__ << ": Filtered short RNA Candidate (< " << minLength << ")\n";
#endif
        continue;
      }
    }

    // remove candidates which don't match the current edge:
    //
    const bool isComplexCand(isComplexSV(readCand));
    if (isComplexCand) {
      if (!readCand.bp1.interval.isIntersect(node1.getInterval())) continue;
      if (!readCand.bp1.interval.isIntersect(node2.getInterval())) continue;
    } else {
      const bool isIntersect(
          (readCand.bp1.interval.isIntersect(node1.getInterval())) &&
          (readCand.bp2.interval.isIntersect(node2.getInterval())));
      const bool isSwapIntersect(
          (readCand.bp1.interval.isIntersect(node2.getInterval())) &&
          (readCand.bp2.interval.isIntersect(node1.getInterval())));
      if (!(isIntersect || isSwapIntersect)) continue;
    }

    // spanning means there's a left|right and left|right breakend pair (in any order) -- note this is not the
    // same as asking if the evidence comes from a read pair. For instance, a CIGAR string can
    // provide a spanning, non-read-pair candidate
    const bool isSpanningCand(isSpanningSV(readCand));

    bool     isMatched(false);
    unsigned svIndex(0);
    for (FatSVCandidate& sv : svs) {
      if (sv.isIntersect(readCand)) {
#if 0
                /// keep candidates formed by semi-mapped reads in separate groups,
                /// these will only be used to augment the evidence of a candidate created with
                /// regular pair evidence -- all purely local candidates will be thrown away.
                ///
                /// the separation starts early (ie. here) because we might not want to use the local-pair
                /// regions... this will take some trial and error
                ///
                const bool isCandLocalOnly(readCand.svEvidenceType == SVEvidenceType::LOCAL_PAIR);
                const bool isSVLocalOnly(sv.bp1.isLocalPairOnly() && sv.bp2.isLocalPairOnly());

                if (isCandLocalOnly == isSVLocalOnly)
#endif
        {
          if (isSpanningCand) {
            // don't store fragment association unless there's a specific hypothesis --
            // if there is no hypothesis (small assembly cases (thus "! isSpanning")), we'll be
            // going back through the bam region during assembly anyway:
            //
            fragment.svLink.emplace_back(svIndex, readCand.svEvidenceType);
          }

          updateEvidenceIndex(fragment, readCand, sv, bamIndex);

          // check evidence distance:
          sv.merge(FatSVCandidate(readCand, bamCount), isExpandSVCandidateSet);

#ifdef DEBUG_SVDATA
          log_os << __FUNCTION__ << ": Added to svIndex: " << svIndex << " match_sv: " << sv << "\n";
#endif

          isMatched = true;
          break;
        }
      }
      svIndex++;
    }

    const bool createNewCandidate(isExpandSVCandidateSet && (!isMatched));
    if (createNewCandidate) {
      const unsigned newSVIndex(svs.size());

#ifdef DEBUG_SVDATA
      log_os << __FUNCTION__ << ": New svIndex: " << newSVIndex << "\n";
#endif

      svs.push_back(FatSVCandidate(readCand, bamCount));
      svs.back().candidateIndex = newSVIndex;

      if (isSpanningCand) {
        // ditto note above, store fragment association only when there's an SV hypothesis:
        fragment.svLink.emplace_back(newSVIndex, readCand.svEvidenceType);
      }
      updateEvidenceIndex(fragment, readCand, svs.back(), bamIndex);
    }
  }
}

void SVFinder::processSequenceFragment(
    const SVLocusNode&              node1,
    const SVLocusNode&              node2,
    const bam_header_info&          bamHeader,
    const reference_contig_segment& node1RefSeq,
    const reference_contig_segment& node2RefSeq,
    const unsigned                  bamIndex,
    const bool                      isExpandSVCandidateSet,
    std::vector<FatSVCandidate>&    svs,
    SVCandidateSetSequenceFragment& fragment,
    SVFinderStats&                  stats)
{
  SVCandidateSetRead* localReadPtr(&(fragment.read1));
  SVCandidateSetRead* remoteReadPtr(&(fragment.read2));
  fragment.svLink.clear();

  if (!localReadPtr->isSet()) {
    std::swap(localReadPtr, remoteReadPtr);
  }

  if (!localReadPtr->isSet()) {
    // this could occur when a supplemental read only is found:
    return;
    //assert(localReadPtr->isSet() && "Neither read in pair is set");
  }

  // sanity check of read pairs
  if (!fragment.checkReadPair()) {
    stats.unmatchedReadPairFilter++;
    return;
  }

  const bam_record* remoteBamRecPtr(remoteReadPtr->isSet() ? &(remoteReadPtr->bamrec) : nullptr);

  const reference_contig_segment& localRef(
      localReadPtr->isSourcedFromGraphEdgeNode1 ? node1RefSeq : node2RefSeq);
  const reference_contig_segment* remoteRefPtr(nullptr);
  if (remoteReadPtr->isSet()) {
    remoteRefPtr = (remoteReadPtr->isSourcedFromGraphEdgeNode1 ? &node1RefSeq : &node2RefSeq);
  }
  _readScanner.getBreakendPair(
      localReadPtr->bamrec, remoteBamRecPtr, bamIndex, bamHeader, localRef, remoteRefPtr, _readCandidates);

  // Label all SVObservations with the sample index:

  // Collapse close spanning sv candidates into complex candidates -- this reflects the fact that the
  // assembler will collapse them anyway, so reduces duplicated work in the assembler;
  for (SVObservation& cand : _readCandidates) {
    if (getSVType(cand) != SV_TYPE::INDEL) continue;
    known_pos_range2   r1(cand.bp1.interval.range);
    known_pos_range2   r2(cand.bp2.interval.range);
    static const pos_t window(30);
    r1.expandBy(window);
    r2.expandBy(window);
    if (!r1.is_range_intersect(r2)) continue;

    // collapse this case:
    cand.bp1.state = SVBreakendState::COMPLEX;
    cand.bp2.state = SVBreakendState::UNKNOWN;
    cand.bp1.interval.range.merge_range(cand.bp2.interval.range);
  }

  // hack split read observations to be symmetrically supported, even though we're only
  // reading in one side:
  for (SVObservation& readCand : _readCandidates) {
    using namespace SVEvidenceType;
    if (readCand.svEvidenceType != SPLIT_ALIGN) continue;

    if (readCand.bp1.lowresEvidence.getVal(SPLIT_ALIGN) == 0) readCand.bp1.lowresEvidence.add(SPLIT_ALIGN);
    if (readCand.bp2.lowresEvidence.getVal(SPLIT_ALIGN) == 0) readCand.bp2.lowresEvidence.add(SPLIT_ALIGN);
  }

#ifdef DEBUG_SVDATA
  log_os << __FUNCTION__ << ": Checking pair: " << fragment << "\n";
  log_os << __FUNCTION__ << ": Translated to candidates:\n";
  for (const SVObservation& cand : _readCandidates) {
    log_os << __FUNCTION__ << ": cand: " << cand << "\n";
  }
#endif
  assignFragmentObservationsToSVCandidates(
      node1, node2, _readCandidates, isExpandSVCandidateSet, fragment, svs, bamIndex);
}

#if 0
static
bool
isLocalEvidence(
    const SVEvidenceType::index_t idx)
{
    using namespace SVEvidenceType;

    switch (idx)
    {
    case CIGAR:
    case SOFTCLIP:
    case SEMIALIGN:
        return true;
    default:
        return false;
    }
}
#endif

static bool isCandidateCountSufficient(const SVCandidate& sv)
{
  static const unsigned           minCandidateComplexCount(2);
  const SVBreakendLowResEvidence& evidence(sv.bp1.lowresEvidence);
  for (unsigned i(0); i < SVEvidenceType::SIZE; ++i) {
    if (SVEvidenceType::isPairType(i)) continue;
    if (evidence.getVal(i) >= minCandidateComplexCount) return true;
  }
  return false;
}

/// Determine if the rate of supporting read observations at a breakpoint is significant
/// relative to a background noise rate
///
/// \param signalReadInfo Vector which has a size equal to the supporting read count, for each supporting
/// read, the vector contains a relative index for the supporting read among all qualifying (mapped) reads
/// from the same input BAM file. This relative index is used to estimate signal density.
///
/// \return true if we reject the null hypothesis that breakpoint signal is noise
static bool isBreakPointSignificant(
    const double alpha, const double noiseRate, std::vector<double>& signalReadInfo)
{
  assert(alpha >= 0.);
  assert((noiseRate >= 0.) && (noiseRate <= 1.));

  const unsigned signalReadCount(signalReadInfo.size());

  // enforce a simple minimum signal count regardless of noiseRate/alpha
  if (signalReadCount < 2) return false;

  assert(signalReadCount >= 1);

  // sort the indices so that by taking the difference we can get an estimate of the
  // number of background reads that occurred 'between' the signal observations:
  std::sort(signalReadInfo.begin(), signalReadInfo.end());
#ifdef DEBUG_SVDATA
  log_os << __FUNCTION__ << ": signalReadInfo_size=" << signalReadInfo.size() << " signalReadInfo={";
  for (const auto rd : signalReadInfo) {
    log_os << rd << ",";
  }
  log_os << "}\n";
#endif

  // The density of signal reads will be variable, so focus on peak density. To do this we pick out a small
  // continuous set of signal observations. The count of observation intervals in this window is
  // 'signalWindowSize' equal to maxSignalWindowSize (or less if the total number of observation
  // intervals is smaller)
  //
  static const unsigned maxSignalWindowSize(4);
  unsigned              signalWindowSize(std::min(maxSignalWindowSize, (signalReadCount - 1)));

  // move a sliding window of size 'signalWindowSize' through the sorted read index list to find the minimum
  // estimated background read count among a continuous signalWindowSize+1 series of signal observations.
  //
  // this greedy minimum approach will tend to err on the side of over-estimating signal rate and calling
  // noisy evidence significant, which is the behavior we want for this kind of conservative noise filter.
  double minWindowBackgroundCount(0);
  for (unsigned signalReadIndex(0); signalReadIndex < (signalReadCount - signalWindowSize);
       ++signalReadIndex) {
    const double windowBackgroundCount(
        signalReadInfo[signalReadIndex + signalWindowSize] - signalReadInfo[signalReadIndex]);
    if ((signalReadIndex == 0) or (minWindowBackgroundCount > windowBackgroundCount)) {
      minWindowBackgroundCount = windowBackgroundCount;
    }
  }
  if (signalWindowSize > minWindowBackgroundCount) signalWindowSize = unsigned(minWindowBackgroundCount);

#ifdef DEBUG_SVDATA
  {
    log_os << __FUNCTION__ << ": noiseRate=" << noiseRate << " signalWindowSize=" << signalWindowSize
           << " backgroundCount=" << minWindowBackgroundCount << " isReject="
           << is_reject_binomial_gte_n_success_exact(
                  alpha, noiseRate, signalWindowSize, static_cast<unsigned>(minWindowBackgroundCount))
           << "\n";
  }
#endif

  return is_reject_binomial_gte_n_success_exact(
      alpha, noiseRate, signalWindowSize, static_cast<unsigned>(minWindowBackgroundCount));
}

/// \brief Test a spanning candidate for minimum supporting evidence level prior
/// to assembly and scoring stages
///
/// Note this test is applied early, and as such it is intended to only filter
/// out cases which are very likely to be noise. This method has an important
/// role in controlling the analysis runtime for FFPE tumor samples, where we
/// might otherwise initiate a very large number of assembly processes for unlikely
/// variant candidates.
static bool isSpanningCandidateSignalSignificant(
    const double noiseRate, const FatSVCandidate& sv, const unsigned bamIndex)
{
  std::vector<double> evidence_bp1;
  std::vector<double> evidence_bp2;
  for (unsigned evidenceTypeIndex(0); evidenceTypeIndex < SVEvidenceType::SIZE; ++evidenceTypeIndex) {
    appendVec(evidence_bp1, sv.bp1EvidenceIndex[evidenceTypeIndex][bamIndex]);
    appendVec(evidence_bp2, sv.bp2EvidenceIndex[evidenceTypeIndex][bamIndex]);
  }

  static const double alpha(0.03);
  const bool          isBp1(isBreakPointSignificant(alpha, noiseRate, evidence_bp1));
  const bool          isBp2(isBreakPointSignificant(alpha, noiseRate, evidence_bp2));

  return (isBp1 || isBp2);
}

static bool isComplexCandidateSignalSignificant(
    const double noiseRate, const FatSVCandidate& sv, const unsigned bamIndex)
{
  std::vector<double> evidence;
  for (unsigned i(0); i < SVEvidenceType::SIZE; ++i) {
    //if (! isLocalEvidence(i)) continue;
    appendVec(evidence, sv.bp1EvidenceIndex[i][bamIndex]);
  }
  static const double alpha(0.005);
  return (isBreakPointSignificant(alpha, noiseRate, evidence));
}

namespace SINGLE_JUNCTION_FILTER {
/// Enumerates different types of filters applied to single sv candidate junctions
enum index_t {
  NONE,
  SEMI_MAPPED,         ///< Candidate is only supported by semi-mapped read pairs
  COMPLEX_LOW_COUNT,   ///< For a complex candidate, absolute evidence read count is too low
  COMPLEX_LOW_SIGNAL,  ///< For a complex candidate, alt allele read rate is not significant wrt expected
                       ///< noise rate
  SPANNING_LOW_SIGNAL  ///< For a spanning candidate, alt allele read rate is not significant wrt expected
                       ///< noise rate
};
}  // namespace SINGLE_JUNCTION_FILTER

static bool isAnySpanningCandidateSignalSignificant(
    const unsigned bamCount, const FatSVCandidate& sv, const std::vector<double>& spanningNoiseRate)
{
  for (unsigned bamIndex(0); bamIndex < bamCount; ++bamIndex) {
    if (isSpanningCandidateSignalSignificant(spanningNoiseRate[bamIndex], sv, bamIndex)) return true;
  }
  return false;
}

static bool isAnyComplexCandidateSignalSignificant(
    const unsigned bamCount, const FatSVCandidate& sv, const std::vector<double>& assemblyNoiseRate)
{
  for (unsigned bamIndex(0); bamIndex < bamCount; ++bamIndex) {
    if (isComplexCandidateSignalSignificant(assemblyNoiseRate[bamIndex], sv, bamIndex)) return true;
  }
  return false;
}

/// Return enumerator value describing the candidate's filtration state, based on
/// information available in a single junction only (as opposed to
/// requiring multi-junction analysis)
///
static SINGLE_JUNCTION_FILTER::index_t isFilterSingleJunctionCandidate(
    const bool                 skipEvidenceSignalFilter,
    const std::vector<double>& spanningNoiseRate,
    const std::vector<double>& assemblyNoiseRate,
    const FatSVCandidate&      sv,
    const unsigned             bamCount)
{
  using namespace SINGLE_JUNCTION_FILTER;

  // don't consider candidates created from only
  // semi-mapped read pairs (ie. one read of the pair is MAPQ0 or MAPQsmall)
  if (sv.bp1.isLocalOnly() && sv.bp2.isLocalOnly()) return SEMI_MAPPED;

  // candidates must have a minimum amount of evidence:
  if (isSpanningSV(sv)) {
    if (!skipEvidenceSignalFilter) {
      if (!isAnySpanningCandidateSignalSignificant(bamCount, sv, spanningNoiseRate))
        return SPANNING_LOW_SIGNAL;
    }
  } else if (isComplexSV(sv)) {
    if (!isCandidateCountSufficient(sv)) return COMPLEX_LOW_COUNT;
    if (!isAnyComplexCandidateSignalSignificant(bamCount, sv, assemblyNoiseRate)) return COMPLEX_LOW_SIGNAL;
  } else {
    assert(false && "Unknown SV candidate type");
  }

  return NONE;
}

/// run early filters on all SV candidates
///
/// some candidates are immediately filtered in this function (so the length of svs will be shorter), other
/// filtration cases require the SV to be marked now and evaluated for filtration once more information is
/// available
static void filterCandidates(
    const bool                   skipEvidenceSignalFilter,
    const std::vector<double>&   spanningNoiseRate,
    const std::vector<double>&   assemblyNoiseRate,
    std::vector<FatSVCandidate>& svs,
    SVFinderStats&               stats,
    const unsigned               bamCount)
{
  unsigned svCount(svs.size());
  unsigned index(0);
  while (index < svCount) {
    using namespace SINGLE_JUNCTION_FILTER;
    const index_t filt(isFilterSingleJunctionCandidate(
        skipEvidenceSignalFilter, spanningNoiseRate, assemblyNoiseRate, svs[index], bamCount));

    bool isFilter(false);
    switch (filt) {
    case SEMI_MAPPED:
      stats.semiMappedFilter++;
      isFilter = true;
      break;
    case COMPLEX_LOW_COUNT:
      stats.ComplexLowCountFilter++;
      isFilter = true;
      break;
    case COMPLEX_LOW_SIGNAL:
      stats.ComplexLowSignalFilter++;
      isFilter = true;
      break;
    case SPANNING_LOW_SIGNAL:
      // this case is marked, but filtering is delayed until a multi-junction evaluation can be done later
      svs[index].isSingleJunctionFilter = true;
      break;
    default:
      break;
    }

    if (isFilter) {
      if ((index + 1) < svCount) svs[index] = svs.back();
      svs.resize(--svCount);
    } else {
      index++;
    }
  }
}

void SVFinder::getCandidatesFromData(
    const SVLocusNode&              node1,
    const SVLocusNode&              node2,
    const bam_header_info&          bamHeader,
    const reference_contig_segment& node1RefSeq,
    const reference_contig_segment& node2RefSeq,
    SVCandidateSetData&             svData,
    std::vector<SVCandidate>&       output_svs,
    SVFinderStats&                  stats)
{
  const unsigned bamCount(_bamStreams.size());

  // track a richer candidates data structure internally, then slice the info down to the
  // regular sv candidate as a last step:
  std::vector<FatSVCandidate> svs;

  for (unsigned bamIndex(0); bamIndex < bamCount; ++bamIndex) {
    SVCandidateSetSequenceFragmentSampleGroup& svDataGroup(svData.getDataGroup(bamIndex));
    for (SVCandidateSetSequenceFragment& fragment : svDataGroup) {
      /// TODO update this test to generalize from read pair to split reads:
      if (!fragment.isAnchored()) continue;

      static const bool isAnchored(true);
      processSequenceFragment(
          node1, node2, bamHeader, node1RefSeq, node2RefSeq, bamIndex, isAnchored, svs, fragment, stats);
    }
  }

  if (_isSomatic) {
    for (unsigned bamIndex(0); bamIndex < bamCount; ++bamIndex) {
      // for somatic calling we're only interested in submapped read processing for the normal sample:
      const bool isTumor(_isAlignmentTumor[bamIndex]);
      if (isTumor) continue;

      SVCandidateSetSequenceFragmentSampleGroup& svDataGroup(svData.getDataGroup(bamIndex));
      for (SVCandidateSetSequenceFragment& pair : svDataGroup) {
        if (pair.isAnchored()) continue;

        static const bool isAnchored(false);
        processSequenceFragment(
            node1, node2, bamHeader, node1RefSeq, node2RefSeq, bamIndex, isAnchored, svs, pair, stats);
      }
    }
  }

#ifdef DEBUG_SVDATA
  {
    log_os << __FUNCTION__ << ": precount: " << svs.size() << "\n";

    unsigned svIndex(0);
    for (SVCandidate& sv : svs) {
      log_os << __FUNCTION__ << ": PRECOUNT: index: " << svIndex << " " << sv;
      svIndex++;
    }
  }
#endif

  consolidateOverlap(bamCount, svData, svs);

#ifdef DEBUG_SVDATA
  {
    log_os << __FUNCTION__ << ": postcount: " << svs.size() << "\n";

    unsigned svIndex(0);
    for (SVCandidate& sv : svs) {
      log_os << __FUNCTION__ << ": POSTCOUNT: index: " << svIndex << " " << sv;
      svIndex++;
    }
  }
#endif

  for (unsigned bamIndex(0); bamIndex < bamCount; ++bamIndex) {
    assert((_spanningNoiseRate[bamIndex] >= 0.) && (_spanningNoiseRate[bamIndex] <= 1.));
    assert((_assemblyNoiseRate[bamIndex] >= 0.) && (_assemblyNoiseRate[bamIndex] <= 1.));
  }

  filterCandidates(_skipEvidenceSignalFilter, _spanningNoiseRate, _assemblyNoiseRate, svs, stats, bamCount);

  std::copy(svs.begin(), svs.end(), std::back_inserter(output_svs));
}

void SVFinder::findCandidateSVImpl(
    const SVLocusSet&         cset,
    const EdgeInfo&           edge,
    SVCandidateSetData&       svData,
    std::vector<SVCandidate>& svs,
    SVFinderStats&            stats)
{
  svData.clear();
  svs.clear();

#ifdef DEBUG_SVDATA
  log_os << "SVDATA: Evaluating edge: " << edge << "\n";
#endif

  // first determine if this is an edge we're going to evaluate
  //
  // edge must be bidirectional at the noise threshold of the locus set:
  if (!isBidirectionalEdge(cset, edge)) {
#ifdef DEBUG_SVDATA
    log_os << "SVDATA: Edge failed min edge count.\n";
#endif
    stats.edgeFilter++;
    return;
  }

  // 1) scan through each region to identify all reads supporting
  // some sort of breakend in the target region, then match up read
  // pairs so that they can easily be accessed from each other
  //
  // 2) iterate through breakend read pairs to estimate the number, type
  // and likely breakend interval regions of SVs corresponding to this edge
  //
  const bam_header_info& bamHeader(cset.getBamHeader());

  const SVLocus& locus(cset.getLocus(edge.locusIndex));

  reference_contig_segment node1RefSeq;
  reference_contig_segment node2RefSeq;
  {
    GenomeInterval searchInterval;
    getNodeRefSeq(bamHeader, locus, edge.nodeIndex1, _referenceFilename, searchInterval, node1RefSeq);
    addSVNodeData(
        bamHeader, locus, edge.nodeIndex1, edge.nodeIndex2, searchInterval, node1RefSeq, true, svData);
  }

  if (edge.nodeIndex1 != edge.nodeIndex2) {
    GenomeInterval searchInterval;
    getNodeRefSeq(bamHeader, locus, edge.nodeIndex2, _referenceFilename, searchInterval, node2RefSeq);
    addSVNodeData(
        bamHeader, locus, edge.nodeIndex2, edge.nodeIndex1, searchInterval, node2RefSeq, false, svData);
  }

  const SVLocusNode& node1(locus.getNode(edge.nodeIndex1));
  const SVLocusNode& node2(locus.getNode(edge.nodeIndex2));
  getCandidatesFromData(node1, node2, bamHeader, node1RefSeq, node2RefSeq, svData, svs, stats);

  //checkResult(svData,svs);
}

void SVFinder::findCandidateSV(
    const SVLocusSet& cset, const EdgeInfo& edge, SVCandidateSetData& svData, std::vector<SVCandidate>& svs)
{
  // time/stats tracking setup:
  const TimeScoper candTime(_edgeTrackerPtr->candidacyTime);
  SVFinderStats    stats;

  findCandidateSVImpl(cset, edge, svData, svs, stats);

  // time/stats tracking finish:
  _edgeStatMan.updateEdgeCandidates(edge, svs.size(), stats);

  if (_isVerbose) {
    log_os << __FUNCTION__
           << ": Low-resolution candidate generation complete. Candidate count: " << svs.size() << "\n";
  }
}
