// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Manta
// Copyright (c) 2013-2014 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/sequencing/licenses/>
//

///
/// \author Chris Saunders
///

#include "SVFinder.hh"

#include "blt_util/bam_streamer.hh"
#include "common/Exceptions.hh"
#include "manta/ReadGroupStatsSet.hh"
#include "manta/SVCandidateUtil.hh"
#include "manta/SVReferenceUtil.hh"

#include "boost/foreach.hpp"

#include <iostream>


// #define DEBUG_SVDATA

#ifdef DEBUG_SVDATA
#include "blt_util/log.hh"
#endif



SVFinder::
SVFinder(const GSCOptions& opt) :
    _scanOpt(opt.scanOpt),
    _isAlignmentTumor(opt.alignFileOpt.isAlignmentTumor),
    _readScanner(_scanOpt,opt.statsFilename,opt.alignFileOpt.alignmentFilename),
    _referenceFilename(opt.referenceFilename),
    _isRNA(opt.isRNA)
{
    // load in set:
    _set.load(opt.graphFilename.c_str(),true);

    _dFilterPtr.reset(new ChromDepthFilterUtil(opt.chromDepthFilename,_scanOpt.maxDepthFactor,_set.header));

    // setup regionless bam_streams:
    // setup all data for main analysis loop:
    BOOST_FOREACH(const std::string& afile, opt.alignFileOpt.alignmentFilename)
    {
        // avoid creating shared_ptr temporaries:
        streamPtr tmp(new bam_streamer(afile.c_str()));
        _bamStreams.push_back(tmp);
    }
}



/// making the dtor explicit and in the cpp allows auto_ptr to work reliably:
SVFinder::
~SVFinder()
{}


// test if read supports an SV on this edge, if so, add to SVData
static
void
addSVNodeRead(
    const std::map<std::string, int32_t>& chromToIndex,
    const SVLocusScanner& scanner,
    const SVLocusNode& localNode,
    const SVLocusNode& remoteNode,
    const bam_record& bamRead,
    const unsigned bamIndex,
    const bool isExpectRepeat,
    const reference_contig_segment& refSeq,
    const bool isNode1,
    SVCandidateSetReadPairSampleGroup& svDataGroup,
    TruthTracker& truthTracker)
{
    using namespace illumina::common;

    if (scanner.isReadFilteredCore(bamRead)) return;
    else if (bamRead.is_unmapped() || bamRead.is_mate_unmapped()) return;

    if (bamRead.map_qual() < scanner.getMinTier2MapQ()) return;

    const bool isNonCompressedAnomalous(scanner.isNonCompressedAnomalous(bamRead,bamIndex));

    bool isLocalAssemblyEvidence(false);
    if (! isNonCompressedAnomalous)
    {
        isLocalAssemblyEvidence = scanner.isLocalAssemblyEvidence(bamRead,refSeq);
    }

    if (! ( isNonCompressedAnomalous || isLocalAssemblyEvidence))
    {
        return; // this read isn't interesting wrt SV discovery
    }

    // finally, check to see if the svDataGroup is full... for now, we allow a very large
    // number of reads to be stored in the hope that we never reach this limit, but just in
    // case we don't want to exhaust memory in centromere pileups, etc...
    //
    // Once svDataGroup is full, we keep scanning but only to find pairs for reads that
    // have already been entered.
    //
    static const unsigned maxDataSize(4000);
    if ((! svDataGroup.isFull()) && (svDataGroup.size() >= maxDataSize))
    {
        svDataGroup.setFull();
    }

    //
    // run an initial screen to make sure at least one candidate from this read matches the regions for this edge:
    //
    typedef std::vector<SVLocus> loci_t;
    loci_t loci;
    scanner.getSVLoci(bamRead, bamIndex, chromToIndex, refSeq, loci,
                      truthTracker);

    BOOST_FOREACH(const SVLocus& locus, loci)
    {
        const unsigned locusSize(locus.size());
        assert((locusSize>=1) && (locusSize<=2));

        unsigned readLocalIndex(0);
        if (locusSize == 2)
        {
            unsigned readRemoteIndex(1);
            if (! locus.getNode(readLocalIndex).isOutCount())
            {
                std::swap(readLocalIndex,readRemoteIndex);
            }

            if (! locus.getNode(readLocalIndex).isOutCount())
            {
                std::ostringstream oss;
                oss << "Unexpected svlocus counts from bam record: " << bamRead << "\n"
                    << "\tlocus: " << locus << "\n";
                BOOST_THROW_EXCEPTION(LogicException(oss.str()));
            }
            if (! locus.getNode(readRemoteIndex).getInterval().isIntersect(remoteNode.getInterval())) continue; //todo should this intersect be checked in swapped orientation?
        }
        else
        {
            if (! locus.getNode(readLocalIndex).getInterval().isIntersect(remoteNode.getInterval())) continue;
        }

        if (! locus.getNode(readLocalIndex).getInterval().isIntersect(localNode.getInterval())) continue; //todo should this intersect be checked in swapped orientation?

        const bool isSubMapped(bamRead.map_qual() < scanner.getMinMapQ());
        svDataGroup.add(bamRead, isExpectRepeat, isNode1, isSubMapped);

        // once any loci has achieved the local/remote overlap criteria, there's no reason to keep scanning loci
        // of the same bam record:
        break;
    }
}



static
void
getNodeRefSeq(
    const bam_header_info& bamHeader,
    const SVLocus& locus,
    const NodeIndexType localNodeIndex,
    const std::string& referenceFilename,
    GenomeInterval& searchInterval,
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
static
void
addReadToDepthEst(
    const bam_record& bamRead,
    const pos_t beginPos,
    std::vector<unsigned>& depth)
{
    const pos_t endPos(beginPos+depth.size());
    const pos_t refStart(bamRead.pos()-1);

    const pos_t readSize(bamRead.read_size());
    for (pos_t readIndex(std::max(0,(beginPos-refStart))); readIndex<readSize; ++readIndex)
    {
        const pos_t refPos(refStart+readIndex);
        if (refPos>=endPos) return;
        const pos_t depthIndex(refPos-beginPos);
        assert(depthIndex>=0);

        depth[depthIndex]++;
    }
}



void
SVFinder::
addSVNodeData(
    const std::map<std::string, int32_t>& chromToIndex,
    const SVLocus& locus,
    const NodeIndexType localNodeIndex,
    const NodeIndexType remoteNodeIndex,
    const GenomeInterval& searchInterval,
    const reference_contig_segment& refSeq,
    const bool isNode1,
    SVCandidateSetData& svData,
    TruthTracker& truthTracker)
{
    // get full search interval:
    const SVLocusNode& localNode(locus.getNode(localNodeIndex));
    const SVLocusNode& remoteNode(locus.getNode(remoteNodeIndex));

    bool isExpectRepeat(svData.setNewSearchInterval(searchInterval));

    /// This is a temporary measure to make the qname collision detection much looser
    /// problems have come up where very large deletions are present in a read, and it is therefore
    /// detected as a repeat in two different regions, even though they are separated by a considerable
    /// distance. Solution is to temporarily turn off collision detection whenever two regions are on
    /// the same chrom (ie. almost always)
    ///
    /// TODO: restore more precise collision detection
    if (! isExpectRepeat) isExpectRepeat = (localNode.getInterval().tid == remoteNode.getInterval().tid);

#ifdef DEBUG_SVDATA
    static const std::string logtag("addSVNodeData: ");
    log_os << logtag << "bp_interval: " << localNode.getInterval()
           << " evidenceInterval: " << localNode.getEvidenceRange()
           << " searchInterval: " << searchInterval
           << " isExpectRepeat: " << isExpectRepeat
           << "\n";
#endif

    const bool isMaxDepth(dFilter().isMaxDepthFilter());
    float maxDepth(0);
    if (isMaxDepth)
    {
        maxDepth = dFilter().maxDepth(searchInterval.tid);
    }
    const pos_t searchBeginPos(searchInterval.range.begin_pos());
    const pos_t searchEndPos(searchInterval.range.end_pos());
    std::vector<unsigned> normalDepthBuffer(searchInterval.range.size(),0);

    bool isFirstTumor(false);

    // iterate through reads, test reads for association and add to svData:
    unsigned bamIndex(0);
    BOOST_FOREACH(streamPtr& bamPtr, _bamStreams)
    {
        const bool isTumor(_isAlignmentTumor[bamIndex]);

        /// assert expected sample order of all normal, then all tumor:
        if (isTumor) isFirstTumor=true;
        assert((! isFirstTumor) || isTumor);

        SVCandidateSetReadPairSampleGroup& svDataGroup(svData.getDataGroup(bamIndex));
        bam_streamer& read_stream(*bamPtr);

        // set bam stream to new search interval:
        read_stream.set_new_region(searchInterval.tid,searchInterval.range.begin_pos(),searchInterval.range.end_pos());

#ifdef DEBUG_SVDATA
        log_os << logtag << "scanning bamIndex: " << bamIndex << "\n";
#endif
        while (read_stream.next())
        {
            const bam_record& bamRead(*(read_stream.get_record_ptr()));

            const pos_t refPos(bamRead.pos()-1);
            if (refPos >= searchEndPos) break;

            if (isMaxDepth)
            {
                if (! isTumor)
                {
                    // depth estimation relies on a simple filtration criteria to stay in sync with the chromosome mean
                    // depth estimates:
                    if (! bamRead.is_unmapped())
                    {
                        addReadToDepthEst(bamRead, searchBeginPos, normalDepthBuffer);
                    }
                }

                assert(refPos<searchEndPos);
                const pos_t depthOffset(refPos - searchBeginPos);
                if ((depthOffset>=0) && (normalDepthBuffer[depthOffset] > maxDepth)) continue;
            }

            // test if read supports an SV on this edge, if so, add to SVData
            addSVNodeRead(
                chromToIndex,_readScanner, localNode, remoteNode,
                bamRead, bamIndex, isExpectRepeat, refSeq, isNode1,
                svDataGroup, truthTracker);
        }
        ++bamIndex;
    }
}



// sanity check the final result
void
SVFinder::
checkResult(
    const SVCandidateSetData& svData,
    const std::vector<SVCandidate>& svs) const
{
    using namespace illumina::common;

    const unsigned svCount(svs.size());
    if (0 == svCount) return;

    // check that the counts totaled up from the data match those in the sv candidates
    std::map<unsigned,unsigned> readCounts;
    std::map<unsigned,unsigned> pairCounts;

    for (unsigned i(0); i<svCount; ++i)
    {
        readCounts[i] = 0;
        pairCounts[i] = 0;
    }

    const unsigned bamCount(_bamStreams.size());
    for (unsigned bamIndex(0); bamIndex < bamCount; ++bamIndex)
    {
        const SVCandidateSetReadPairSampleGroup& svDataGroup(svData.getDataGroup(bamIndex));
        BOOST_FOREACH(const SVCandidateSetReadPair& pair, svDataGroup)
        {
            BOOST_FOREACH(const SVPairAssociation& sva, pair.svLink)
            {
                if (sva.index>=svCount)
                {
                    std::ostringstream oss;
                    oss << "Searching for SVIndex: " << sva.index << " with svSize: " << svCount << "\n";
                    BOOST_THROW_EXCEPTION(LogicException(oss.str()));
                }

                if (SVEvidenceType::isPairType(sva.evtype))
                {
                    if (pair.read1.isSet()) readCounts[sva.index]++;
                    if (pair.read2.isSet()) readCounts[sva.index]++;
                    if (pair.read1.isSet() && pair.read2.isSet()) pairCounts[sva.index] += 2;
                }
            }
        }
    }

    for (unsigned svIndex(0); svIndex<svCount; ++svIndex)
    {
        const unsigned svObsReadCount(svs[svIndex].bp1.getLocalPairCount() + svs[svIndex].bp2.getLocalPairCount());
        const unsigned svObsPairCount(svs[svIndex].bp1.getPairCount() + svs[svIndex].bp2.getPairCount());
        assert(svs[svIndex].bp1.getPairCount() == svs[svIndex].bp2.getPairCount());

        const unsigned dataObsReadCount(readCounts[svIndex]);
        const unsigned dataObsPairCount(pairCounts[svIndex]);

        bool isCountException(false);
        if      (svObsReadCount != dataObsReadCount) isCountException=true;
        else if (svObsPairCount != dataObsPairCount) isCountException=true;

        if (isCountException)
        {
            std::ostringstream oss;
            oss << "Unexpected difference in sv and data read counts.\n"
                << "\tSVreadCount: " << svObsReadCount << " DataReadCount: " << dataObsReadCount << "\n"
                << "\tSVpaircount: " << svObsPairCount << " DataPaircount: " << dataObsPairCount << "\n"
                << "\tsvIndex: " << svIndex << " SV: " << svs[svIndex];
            BOOST_THROW_EXCEPTION(LogicException(oss.str()));

        }
    }
}



typedef std::map<unsigned,unsigned> movemap_t;



/// local convenience struct, if only I had closures instead... :<
struct svCandDeleter
{
    svCandDeleter(
        std::vector<SVCandidate>& svs,
        movemap_t& moveSVIndex) :
        _shift(0),
        _isLastIndex(false),
        _lastIndex(0),
        _svs(svs),
        _moveSVIndex(moveSVIndex)
    {}

    void
    deleteIndex(
        const unsigned index)
    {
        assert(index <= _svs.size());

        if (_isLastIndex)
        {
            for (unsigned i(_lastIndex+1); i<index; ++i)
            {
                assert(_shift>0);
                assert(i>=_shift);

                _svs[(i-_shift)] = _svs[i];
                // moveSVIndex has already been set for deleted indices, this sets
                // the move for non-deleted positions:
                _moveSVIndex[i] = (i-_shift);
            }
        }
        _lastIndex=index;
        _isLastIndex=true;
        _shift++;
    }

private:
    unsigned _shift;
    bool _isLastIndex;
    unsigned _lastIndex;
    std::vector<SVCandidate>& _svs;
    movemap_t& _moveSVIndex;
};



// check whether any svs have grown to intersect each other
//
// this is also part of the temp hygen hack, so just make this minimally work:
//
static
void
consolidateOverlap(
    const unsigned bamCount,
    SVCandidateSetData& svData,
    std::vector<SVCandidate>& svs)
{
#ifdef DEBUG_SVDATA
    static const std::string logtag("consolidateOverlap: ");
#endif

    movemap_t moveSVIndex;
    std::set<unsigned> deletedSVIndex;

    std::vector<unsigned> innerIndexShift;

    const unsigned svCount(svs.size());
    for (unsigned outerIndex(1); outerIndex<svCount; ++outerIndex)
    {
        const unsigned prevInnerIndexShift( (outerIndex<=1) ? 0 : innerIndexShift[outerIndex-2]);
        innerIndexShift.push_back(prevInnerIndexShift + deletedSVIndex.count(outerIndex-1));
        for (unsigned innerIndex(0); innerIndex<outerIndex; ++innerIndex)
        {
            if (deletedSVIndex.count(innerIndex)) continue;

            if (svs[innerIndex].isIntersect(svs[outerIndex]))
            {
#ifdef DEBUG_SVDATA
                log_os << logtag << "Merging outer:inner: " << outerIndex << " " << innerIndex << "\n";
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

    if (! deletedSVIndex.empty())
    {
#ifdef DEBUG_SVDATA
        BOOST_FOREACH(const unsigned index, deletedSVIndex)
        {
            log_os << logtag << "deleted index: " << index << "\n";
        }
#endif

        {
            svCandDeleter svDeleter(svs,moveSVIndex);

            BOOST_FOREACH(const unsigned index, deletedSVIndex)
            {
                svDeleter.deleteIndex(index);
            }
            svDeleter.deleteIndex(svCount);
        }

        svs.resize(svs.size()-deletedSVIndex.size());

        // fix indices:
        for (unsigned i(0); i<svs.size(); ++i)
        {
            svs[i].candidateIndex = i;
        }
    }

    if (! moveSVIndex.empty())
    {
#ifdef DEBUG_SVDATA
        BOOST_FOREACH(const movemap_t::value_type& val, moveSVIndex)
        {
            log_os << logtag << "Movemap from: " << val.first << " to: " << val.second << "\n";
        }
#endif

        for (unsigned bamIndex(0); bamIndex < bamCount; ++bamIndex)
        {
            SVCandidateSetReadPairSampleGroup& svDataGroup(svData.getDataGroup(bamIndex));
            BOOST_FOREACH(SVCandidateSetReadPair& pair, svDataGroup)
            {
                BOOST_FOREACH(SVPairAssociation& sva, pair.svLink)
                {
                    if (moveSVIndex.count(sva.index))
                    {
                        sva.index = moveSVIndex[sva.index];
                    }
                }
            }
        }
    }
}


/// readCandidates are the set of hypotheses generated by individual read pair --
/// this is the read pair which we seek to assign to one of the identified SVs (in svs)
/// or we push the candidate into svs to start a new candidate associated with this edge
///
/// this is meant as only a temporary form of hypothesis generation, in the current system
/// we do at least delineate alternative candidates by strand and region overlap, but over
/// the longer term we should be able to delineate cluster by a clustering of possible
/// breakend locations.
///
/// isExpandSVCandidateSet if false, don't add new SVs or expand existing SVs
///
void
SVFinder::
assignPairObservationsToSVCandidates(
    const SVLocusNode& node1,
    const SVLocusNode& node2,
    const std::vector<SVObservation>& readCandidates,
    const bool isExpandSVCandidateSet,
    SVCandidateSetReadPair& pair,
    std::vector<SVCandidate>& svs)
{
#ifdef DEBUG_SVDATA
    static const std::string logtag("assignPairObservationsToSVCandidates: ");
#endif

    // we anticipate so few svs from the POC method, that there's no indexing on them
    // OST 26/09/2013: Be careful when re-arranging or rewriting the code below, under g++ 4.1.2
    // this can lead to an infinite loop.
    BOOST_FOREACH(const SVObservation& readCand, readCandidates)
    {
#ifdef DEBUG_SVDATA
        log_os << logtag << "Starting assignment for read cand: " << readCand << "\n";
#endif
        if (_isRNA)
        {
            if ((! SVBreakendState::isSameOrientation(readCand.bp1.state,readCand.bp2.state) ||
                 (! SVBreakendState::isSimpleBreakend(readCand.bp1.state)) ||
                 (! SVBreakendState::isSimpleBreakend(readCand.bp2.state))) &&
                isSVBelowMinSize(readCand, _scanOpt.minRNALength))
            {
#ifdef DEBUG_SVDATA
                log_os << logtag << "Filtered short RNA Candidate\n";
#endif
                continue;
            }
        }

        // remove candidates which don't match the current edge:
        //
        if (isComplexSV(readCand))
        {
            if (! readCand.bp1.interval.isIntersect(node1.getInterval())) continue;
            if (! readCand.bp1.interval.isIntersect(node2.getInterval())) continue;
        }
        else
        {
            const bool isIntersect((readCand.bp1.interval.isIntersect(node1.getInterval())) &&
                                   (readCand.bp2.interval.isIntersect(node2.getInterval())));
            const bool isSwapIntersect((readCand.bp1.interval.isIntersect(node2.getInterval())) &&
                                       (readCand.bp2.interval.isIntersect(node1.getInterval())));
            if (! (isIntersect || isSwapIntersect)) continue;
        }

        /// spanning means there's a left|right and left|right breakend pair (in any order) -- note this is not the
        /// same as asking if the evidence comes from a read pair. For instance, a CIGAR string can
        /// provide a spanning, non-read-pair candidate
        const bool isSpanning(isSpanningSV(readCand));

        bool isMatched(false);
        unsigned svIndex(0);
        BOOST_FOREACH(SVCandidate& sv, svs)
        {
            if (sv.isIntersect(readCand))
            {
#if 0
                /// keep candidates formed by semi-mapped reads in separate groups,
                /// these will only be used to augment the evidence of a candidate created with
                /// regular pair evidence -- all purely local candidates will be thrown away.
                ///
                /// the separation starts early (ie. here) because we might not want to use the local-pair
                /// regions... this will take some trial and error
                ///
                const bool isCandLocalOnly(readCand.evtype == SVEvidenceType::LOCAL_PAIR);
                const bool isSVLocalOnly(sv.bp1.isLocalPairOnly() && sv.bp2.isLocalPairOnly());

                if (isCandLocalOnly == isSVLocalOnly)
#endif
                {

#ifdef DEBUG_SVDATA
                    log_os << logtag << "Adding to svIndex: " << svIndex << " match_sv: " << sv << "\n";
#endif
                    if (isSpanning)
                    {
                        // don't store read pairs unless there's a specific hypothesis --
                        // if there is no hypothesis (small assembly cases (thus "! isSpanning")), we'll be
                        // going back through the bam region during assembly anyway:
                        //
                        pair.svLink.push_back(SVPairAssociation(svIndex,readCand.evtype));
                    }

                    if (isExpandSVCandidateSet)
                    {
                        sv.merge(readCand);
                    }
                    else
                    {
                        // add submapped read pair evidence -- but don't allow these reads to expand the current candidate size:
                        sv.evidenceMerge(readCand);
                    }

                    isMatched=true;
                    break;
                }
            }
            svIndex++;
        }

        const bool createNewCandidate(isExpandSVCandidateSet && (! isMatched));
        if (createNewCandidate)
        {
            const unsigned newSVIndex(svs.size());

#ifdef DEBUG_SVDATA
            log_os << logtag << "New svIndex: " << newSVIndex << "\n";
#endif
            if (isSpanning)
            {
                // ditto note above, store read pairs only when there's an SV hypothesis:
                pair.svLink.push_back(SVPairAssociation(newSVIndex,readCand.evtype));
            }
            svs.push_back(readCand);
            svs.back().candidateIndex = newSVIndex;
        }
    }
}



/// we either process the pair to discover new SVs and expand existing SVs,
/// or we go through and add pairs to existing SVs without expansion
///
void
SVFinder::
processReadPair(
    const SVLocusNode& node1,
    const SVLocusNode& node2,
    const std::map<std::string, int32_t>& chromToIndex,
    const reference_contig_segment& refSeq1,
    const reference_contig_segment& refSeq2,
    const unsigned bamIndex,
    const bool isExpandSVCandidateSet,
    std::vector<SVCandidate>& svs,
    TruthTracker& truthTracker,
    SVCandidateSetReadPair& pair)
{
     SVCandidateSetRead* localReadPtr(&(pair.read1));
     SVCandidateSetRead* remoteReadPtr(&(pair.read2));
     pair.svLink.clear();

     if (! localReadPtr->isSet())
     {
         std::swap(localReadPtr,remoteReadPtr);
     }
     assert(localReadPtr->isSet() && "Neither read in pair is set");

     const bam_record* remoteBamRecPtr( remoteReadPtr->isSet() ? &(remoteReadPtr->bamrec) : NULL);

     const reference_contig_segment& localRef( localReadPtr->isNode1 ? refSeq1 : refSeq2 );
     const reference_contig_segment* remoteRefPtr(NULL);
     if (remoteReadPtr->isSet())
     {
         remoteRefPtr = (remoteReadPtr->isNode1 ?  &refSeq1 : &refSeq2 );
     }
     _readScanner.getBreakendPair(localReadPtr->bamrec, remoteBamRecPtr,
                                  bamIndex, chromToIndex, localRef,
                                  remoteRefPtr, _readCandidates,
                                  truthTracker);

#ifdef DEBUG_SVDATA
     log_os << __FUNCTION__ << ": Checking pair: " << pair << "\n";
     log_os << __FUNCTION__ << ": Translated to candidates:\n";
     BOOST_FOREACH(const SVObservation& cand, readCandidates)
     {
         log_os << __FUNCTION__ << ": cand: " << cand << "\n";
     }
#endif
     assignPairObservationsToSVCandidates(node1, node2, _readCandidates, isExpandSVCandidateSet, pair, svs);
}


void
SVFinder::
getCandidatesFromData(
    const SVLocusNode& node1,
    const SVLocusNode& node2,
    const std::map<std::string, int32_t>& chromToIndex,
    const reference_contig_segment& refSeq1,
    const reference_contig_segment& refSeq2,
    SVCandidateSetData& svData,
    std::vector<SVCandidate>& svs,
    TruthTracker& truthTracker)
{
    const unsigned bamCount(_bamStreams.size());

    for (unsigned bamIndex(0); bamIndex<bamCount; ++bamIndex)
    {
        SVCandidateSetReadPairSampleGroup& svDataGroup(svData.getDataGroup(bamIndex));
        BOOST_FOREACH(SVCandidateSetReadPair& pair, svDataGroup)
        {
            if (! pair.isAnchored()) continue;

            static const bool isAnchored(true);
            processReadPair(
                node1, node2, chromToIndex, refSeq1, refSeq2, bamIndex, isAnchored,
                svs, truthTracker, pair);
        }
    }

    for (unsigned bamIndex(0); bamIndex<bamCount; ++bamIndex)
    {
        SVCandidateSetReadPairSampleGroup& svDataGroup(svData.getDataGroup(bamIndex));
        BOOST_FOREACH(SVCandidateSetReadPair& pair, svDataGroup)
        {
            if (pair.isAnchored()) continue;

            static const bool isAnchored(false);
            processReadPair(
                node1, node2, chromToIndex, refSeq1, refSeq2, bamIndex, isAnchored,
                svs, truthTracker, pair);
        }
    }

#ifdef DEBUG_SVDATA
    {
        log_os << __FUNCTION__ << ": precount: " << svs.size() << "\n";

        unsigned svIndex(0);
        BOOST_FOREACH(SVCandidate& sv, svs)
        {
            log_os << __FUNCTION__ << ": PRECOUNT: index: " << svIndex << " " << sv;
            svIndex++;
        }
    }
#endif

    consolidateOverlap(bamCount,svData,svs);

#ifdef DEBUG_SVDATA
    {
        log_os << __FUNCTION__ << ": postcount: " << svs.size() << "\n";

        unsigned svIndex(0);
        BOOST_FOREACH(SVCandidate& sv, svs)
        {
            log_os << __FUNCTION__ << ": POSTCOUNT: index: " << svIndex << " " << sv;
            svIndex++;
        }
    }
#endif
}



void
SVFinder::
findCandidateSV(
    const std::map<std::string, int32_t>& chromToIndex,
    const EdgeInfo& edge,
    SVCandidateSetData& svData,
    std::vector<SVCandidate>& svs,
    TruthTracker& truthTracker)
{
    svData.clear();
    svs.clear();

    const SVLocusSet& set(getSet());
    const unsigned minEdgeCount(set.getMinMergeEdgeCount());

#ifdef DEBUG_SVDATA
    log_os << "SVDATA: Evaluating edge: " << edge << "\n";
#endif

    // first determine if this is an edge we're going to evaluate
    //
    // edge must be bidirectional at the noise threshold of the locus set:
    const SVLocus& locus(set.getLocus(edge.locusIndex));

    if ((locus.getEdge(edge.nodeIndex1,edge.nodeIndex2).getCount() < minEdgeCount) ||
        (locus.getEdge(edge.nodeIndex2,edge.nodeIndex1).getCount() < minEdgeCount))
    {
#ifdef DEBUG_SVDATA
        log_os << "SVDATA: Edge failed min edge count.\n";
#endif
        return;
    }

    //
    // 1) scan through each region to identify all reads supporting
    // some sort of breakend in the target region, then match up read
    // pairs so that they can easily be accessed from each other
    //
    // 2) iterate through breakend read pairs to estimate the number, type
    // and likely breakend interval regions of SVs corresponding to this edge
    //
    const bam_header_info& bamHeader(set.header);

    reference_contig_segment refSeq1;
    reference_contig_segment refSeq2;
    {
        GenomeInterval searchInterval;
        getNodeRefSeq(bamHeader, locus, edge.nodeIndex1, _referenceFilename, searchInterval, refSeq1);
        addSVNodeData(chromToIndex, locus, edge.nodeIndex1, edge.nodeIndex2,
                      searchInterval, refSeq1, true, svData, truthTracker);
    }

    if (edge.nodeIndex1 != edge.nodeIndex2)
    {
        GenomeInterval searchInterval;
        getNodeRefSeq(bamHeader, locus, edge.nodeIndex2, _referenceFilename, searchInterval, refSeq2);
        addSVNodeData(chromToIndex, locus, edge.nodeIndex2, edge.nodeIndex1,
                      searchInterval, refSeq2, false, svData, truthTracker);
    }

    const SVLocusNode& node1(locus.getNode(edge.nodeIndex1));
    const SVLocusNode& node2(locus.getNode(edge.nodeIndex2));
    getCandidatesFromData(node1,node2,chromToIndex, refSeq1, refSeq2, svData,
                          svs, truthTracker);

    //checkResult(svData,svs);
}

