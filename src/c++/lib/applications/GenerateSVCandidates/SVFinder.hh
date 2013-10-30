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
// <https://github.com/sequencing/licenses/>
//

///
/// \author Chris Saunders
///

#pragma once

#include "GSCOptions.hh"
#include "svgraph/EdgeInfo.hh"

#include "blt_util/bam_streamer.hh"
#include "manta/SVCandidate.hh"
#include "manta/SVCandidateSetData.hh"
#include "manta/SVLocusScanner.hh"
#include "svgraph/SVLocusSet.hh"

#include "boost/shared_ptr.hpp"

#include <vector>


struct SVFinder
{

    SVFinder(const GSCOptions& opt);

    const SVLocusSet&
    getSet() const
    {
        return _set;
    }

    void
    findCandidateSV(
        const std::map<std::string, int32_t>& chromToIndex,
        const EdgeInfo& edge,
        SVCandidateSetData& svData,
        std::vector<SVCandidate>& svs);

    void
    checkResult(
        const SVCandidateSetData& svData,
        const std::vector<SVCandidate>& svs) const;

private:

    void
    addSVNodeData(
        const std::map<std::string, int32_t>& chromToIndex,
        const SVLocus& locus,
        const NodeIndexType node1,
        const NodeIndexType node2,
        SVCandidateSetData& svData);


    void
    getCandidatesFromData(
        const std::map<std::string, int32_t>& chromToIndex,
        SVCandidateSetData& svData,
        std::vector<SVCandidate>& svs);

    const ReadScannerOptions _scanOpt;
    SVLocusSet _set;
    SVLocusScanner _readScanner;

    typedef boost::shared_ptr<bam_streamer> streamPtr;
    std::vector<streamPtr> _bamStreams;
};
