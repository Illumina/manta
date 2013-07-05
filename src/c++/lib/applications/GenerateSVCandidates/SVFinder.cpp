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
#include "manta/SVCandidateCachedData.hh"

#include "boost/foreach.hpp"

#include <iostream>



SVFinder::
SVFinder(const GSCOptions& opt)
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
    SVCandidateCachedData svData;
}



