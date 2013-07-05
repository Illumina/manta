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

#include "blt_util/input_stream_handler.hh"
#include "blt_util/bam_streamer.hh"
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
        SVCandidateData& svData)
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

    svData.add(read,bamIndex);
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

    // set all bam streams to new search interval:
    BOOST_FOREACH(streamPtr& bamPtr, _bamStreams)
    {
        bamPtr->set_new_region(searchInterval.tid,searchInterval.range.begin_pos(),searchInterval.range.end_pos());
    }

    // iterate through reads, test reads for association and add to svData:
    unsigned bamIndex(0);
    BOOST_FOREACH(streamPtr& bamPtr, _bamStreams)
    {
        bam_streamer& read_stream(*bamPtr);
        while (read_stream.next())
        {
            const bam_record& read(*(read_stream.get_record_ptr()));

            // test if read supports an SV on this edge, if so, add to SVData
            addSVNodeRead(_readScanner,localNode,remoteNode,read, bamIndex,svData);
        }
        bamIndex++;
    }
}



void
SVFinder::
findSVCandidates(
        const EdgeInfo& edge,
        std::vector<SVCandidate>& svs)
{
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

    SVCandidateData svData;

    addSVNodeData(locus,edge.nodeIndex1,edge.nodeIndex2,svData);
    addSVNodeData(locus,edge.nodeIndex1,edge.nodeIndex2,svData);
}



