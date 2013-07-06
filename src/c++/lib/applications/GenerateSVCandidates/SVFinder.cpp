// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Manta
// Copyright (c) 2013 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/downloads/sequencing/licenses/>.
//

///
/// \author Chris Saunders
///

#include "SVFinder.hh"

#include "blt_util/bam_streamer.hh"
#include "blt_util/input_stream_handler.hh"
#include "blt_util/log.hh"
#include "manta/ReadGroupStatsSet.hh"

#include "boost/foreach.hpp"

#include <iostream>



SVFinder::
SVFinder(const GSCOptions& opt) :
    _readScanner(opt.scanOpt,opt.statsFilename,opt.alignmentFilename)
{
    // load in set:
    _set.load(opt.graphFilename.c_str());

    // setup regionless bam_streams:
    // setup all data for main analysis loop:
    BOOST_FOREACH(const std::string& afile, opt.alignmentFilename)
    {
        // avoid creating shared_ptr temporaries:
        streamPtr tmp(new bam_streamer(afile.c_str()));
        _bamStreams.push_back(tmp);
    }
}



// test if read supports an SV on this edge, if so, add to SVData
static
void
addSVNodeRead(
        const SVLocusScanner& scanner,
        const SVLocusNode& localNode,
        const SVLocusNode& remoteNode,
        const bam_record& read,
        const unsigned bamIndex,
        SVCandidateDataGroup& svDataGroup)
{
    if(scanner.isReadFiltered(read)) return;

    SVLocus locus;
    scanner.getChimericSVLocus(read,bamIndex,locus);
    const SVLocus& clocus(locus);

    if(clocus.empty()) return;
    if(2 != clocus.size()) return;

    unsigned readLocalIndex(0);
    unsigned readRemoteIndex(1);
    if(0 == clocus.getNode(readLocalIndex).count)
    {
        std::swap(readLocalIndex,readRemoteIndex);
    }

    if(! clocus.getNode(readLocalIndex).interval.isIntersect(localNode.interval)) return;
    if(! clocus.getNode(readRemoteIndex).interval.isIntersect(remoteNode.interval)) return;

    svDataGroup.add(read);
}



void
SVFinder::
addSVNodeData(
        const SVLocus& locus,
        const NodeIndexType localNodeIndex,
        const NodeIndexType remoteNodeIndex,
        SVCandidateData& svData)
{
    // get full search interval:
    const SVLocusNode& localNode(locus.getNode(localNodeIndex));
    const SVLocusNode& remoteNode(locus.getNode(remoteNodeIndex));
    GenomeInterval searchInterval(localNode.interval);

    searchInterval.range.merge_range(localNode.evidenceRange);

#ifdef DEBUG_SVDATA
    log_os << "addSVNodeData: bp_interval: " << localNode.interval
           << " evidnece_interval: " << localNode.evidenceRange
           << " search_interval: " << searchInterval
           << "\n";
#endif

    // iterate through reads, test reads for association and add to svData:
    unsigned bamIndex(0);
    BOOST_FOREACH(streamPtr& bamPtr, _bamStreams)
    {
        SVCandidateDataGroup& svDataGroup(svData.getDataGroup(bamIndex));
        bam_streamer& read_stream(*bamPtr);

        // set bam stream to new search interval:
        read_stream.set_new_region(searchInterval.tid,searchInterval.range.begin_pos(),searchInterval.range.end_pos());

#ifdef DEBUG_SVDATA
    log_os << "addSVNodeData: scanning bamIndex: " << bamIndex << "\n";
#endif
        while (read_stream.next())
        {
            const bam_record& bamRead(*(read_stream.get_record_ptr()));

            // test if read supports an SV on this edge, if so, add to SVData
            addSVNodeRead(_readScanner,localNode,remoteNode,bamRead, bamIndex,svDataGroup);
        }
        bamIndex++;
    }
}



// sanity check the final result
void
SVFinder::
checkResult(
        const SVCandidateData& svData,
        const std::vector<SVCandidate>& svs) const
{
    // check that the counts totalled up from the data match those in the sv candidates
    std::map<unsigned,unsigned> counts;
    const unsigned svCount(svs.size());
    for(unsigned i(0);i<svCount;++i)
    {
        counts[i] = 0;
    }

    const unsigned bamCount(_bamStreams.size());
    for(unsigned bamIndex(0); bamIndex < bamCount; ++bamIndex)
    {
        const SVCandidateDataGroup& svDataGroup(svData.getDataGroup(bamIndex));
        BOOST_FOREACH(const SVCandidateReadPair& pair, svDataGroup)
        {
            assert(pair.svIndex<svCount);
            if(pair.read1.isSet()) counts[pair.svIndex]++;
            if(pair.read2.isSet()) counts[pair.svIndex]++;
        }
    }

    for(unsigned svIndex(0); svIndex<svCount; ++svIndex)
    {
        const unsigned svObsCount(svs[svIndex].bp1.readCount + svs[svIndex].bp2.readCount);
        const unsigned dataObsCount(counts[svIndex]);
        assert(svObsCount == dataObsCount);
    }
}


void
SVFinder::
getCandidatesFromData(
        SVCandidateData& svData,
        std::vector<SVCandidate>& svs)
{
    SVCandidate cand;

    const unsigned bamCount(_bamStreams.size());
    for(unsigned bamIndex(0); bamIndex < bamCount; ++bamIndex)
    {
        SVCandidateDataGroup& svDataGroup(svData.getDataGroup(bamIndex));
        BOOST_FOREACH(SVCandidateReadPair& pair, svDataGroup)
        {
            SVCandidateRead* localReadPtr(&(pair.read1));
            SVCandidateRead* remoteReadPtr(&(pair.read2));
            if(! localReadPtr->isSet())
            {
                std::swap(localReadPtr,remoteReadPtr);
            }
            const bam_record* remoteBamRecPtr( remoteReadPtr->isSet() ? &(remoteReadPtr->bamrec) : NULL);

            cand.clear();
            _readScanner.getBreakendPair(localReadPtr->bamrec, remoteBamRecPtr, bamIndex, cand.bp1, cand.bp2);

#ifdef DEBUG_SVDATA
            log_os << "Checking pair: " << pair << "\n";
            log_os << "Translated to cand: " << cand << "\n";
#endif

            bool isSVFound(false);
            unsigned svIndex(0);

            // temporary hack hypoth gen method assumes that only one SV exists for each overlapping breakpoint range with
            // the same orientation:
            //
            // we anticipate so few svs from the POC method, that there's no indexing on them
            BOOST_FOREACH(SVCandidate& sv, svs)
            {
                if(sv.isIntersect(cand))
                {
#ifdef DEBUG_SVDATA
                    log_os << "Adding to svIndex: " << svIndex << "\n";
#endif
                    sv.merge(cand);
                    pair.svIndex = svIndex;
                    isSVFound=true;
                    break;
                }
                svIndex++;
            }

            if(! isSVFound)
            {
#ifdef DEBUG_SVDATA
                    log_os << "New svIndex: " << svs.size() << "\n";
#endif
                pair.svIndex = svs.size();
                svs.push_back(cand);
            }
        }
    }

#ifdef DEBUG_SVDATA
    log_os << "findSVCandidates: precount: " << svs.size() << "\n";
    BOOST_FOREACH(SVCandidate& sv, svs)
    {
        log_os << "\tPRECOUNT: " << sv << "\n";
    }
#endif


    // finally check whether any svs have grown to intersect each other
    //
    // this is also part of the temp hygen hack, so just make it function:
    //
    std::map<unsigned,unsigned> moveSVIndex;
    std::set<unsigned> deletedSVIndex;
    const unsigned svCount(svs.size());
    for(unsigned outerIndex(0); outerIndex<svCount; ++outerIndex)
    {
        const unsigned routerIndex(svCount-(outerIndex+1));
        for(unsigned innerIndex(0); innerIndex<routerIndex; ++innerIndex)
        {
            if(svs[innerIndex].isIntersect(svs[routerIndex]))
            {
#ifdef DEBUG_SVDATA
                log_os << "Merging outer:inner: " << routerIndex << " " << innerIndex << "\n";
#endif
                svs[innerIndex].merge(svs[routerIndex]);
                moveSVIndex[routerIndex] = innerIndex;
                deletedSVIndex.insert(routerIndex);
                break;
            }
        }
    }

    if(! deletedSVIndex.empty())
    {
        {
            unsigned shift(0);
            bool isLastIndex(false);
            unsigned lastIndex(0);
            BOOST_FOREACH(const unsigned index, deletedSVIndex)
            {
                shift++;
                if(isLastIndex)
                {
                    for(unsigned i(lastIndex+1);i<index;++i)
                    {
                        assert(shift>0);
                        assert(i>=shift);
                        moveSVIndex[i] = (i-shift);
                        svs[(i-shift)] = svs[i];
                    }
                }
                lastIndex=index;
                isLastIndex=true;
            }
            if(isLastIndex)
            {
                for(unsigned i(lastIndex+1);i<svCount;++i)
                {
                    assert(shift>0);
                    assert(i>=shift);
                    moveSVIndex[i] = (i-shift);
                    svs[(i-shift)] = svs[i];
                }
            }
        }

        svs.resize(svs.size()-deletedSVIndex.size());
    }

    if(! moveSVIndex.empty())
    {
        for(unsigned bamIndex(0); bamIndex < bamCount; ++bamIndex)
        {
            SVCandidateDataGroup& svDataGroup(svData.getDataGroup(bamIndex));
            BOOST_FOREACH(SVCandidateReadPair& pair, svDataGroup)
            {
                if(moveSVIndex.count(pair.svIndex))
                {
                    pair.svIndex = moveSVIndex[pair.svIndex];
                }
            }
        }
    }

#ifdef DEBUG_SVDATA
    log_os << "findSVCandidates: postcount: " << svs.size() << "\n";
    BOOST_FOREACH(SVCandidate& sv, svs)
    {
        log_os << "\tPOSTCOUNT: " << sv << "\n";
    }
#endif
}



void
SVFinder::
findSVCandidates(
        const EdgeInfo& edge,
        SVCandidateData& svData,
        std::vector<SVCandidate>& svs)
{
    svData.clear();
    svs.clear();

    const SVLocusSet& set(getSet());
    const unsigned minEdgeCount(set.getMinMergeEdgeCount());

    // first determine if this is an edge we're going to evaluate
    //
    // edge must be bidirectional at the noise threshold of the locus set:
    const SVLocus& locus(set.getLocus(edge.locusIndex));

    if((locus.getEdge(edge.nodeIndex1,edge.nodeIndex2).count <= minEdgeCount) ||
       (locus.getEdge(edge.nodeIndex2,edge.nodeIndex1).count <= minEdgeCount))
    {
        return;
    }

    // start gathering evidence required for hypothesis generation,
    //
    // first step is to scan through each region, and identify all reads supporting
    // some sort of breakend in the target region, then match up pairs so that they
    // can easily be accessed from each other
    //

    // steps:
    // iterate through regions -- for each region walk from evidence range to breakpoint range picking up all reads associated with the breakend
    // make a common code to determine read breakend association shared with Estimation step
    // pair reads in data structure
    // determine the number of breakends from simple orientation and intersecting joint region logic
    // assign data to each breakend candidates
    // come up with an ultra-simple model-free scoring rule: >10 obs = Q60,k else Q0
    //

    addSVNodeData(locus,edge.nodeIndex1,edge.nodeIndex2,svData);
    addSVNodeData(locus,edge.nodeIndex2,edge.nodeIndex1,svData);

    getCandidatesFromData(svData,svs);

#ifdef DEBUG_SVDATA
    checkResult(svData,svs);
#endif
}



