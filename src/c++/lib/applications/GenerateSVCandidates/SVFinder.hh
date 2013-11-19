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

#include "blt_util/bam_streamer.hh"
#include "manta/SVCandidate.hh"
#include "manta/SVCandidateSetData.hh"
#include "manta/SVLocusScanner.hh"
#include "svgraph/EdgeInfo.hh"
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
        std::vector<SVCandidate>& svs,
        TruthTracker& truthTracker);

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
        const GenomeInterval& searchInterval,
        const reference_contig_segment& refSeq,
        const bool isNode1,
        SVCandidateSetData& svData,
        TruthTracker& truthTracker);


    void
    getCandidatesFromData(
        const SVLocusNode& node1,
        const SVLocusNode& node2,
        const std::map<std::string, int32_t>& chromToIndex,
        const reference_contig_segment& refSeq1,
        const reference_contig_segment& refSeq2,
        SVCandidateSetData& svData,
        std::vector<SVCandidate>& svs,
        TruthTracker& truthTracker);

    const ReadScannerOptions _scanOpt;
    SVLocusSet _set;
    SVLocusScanner _readScanner;

    const std::string _referenceFilename;

    typedef boost::shared_ptr<bam_streamer> streamPtr;
    std::vector<streamPtr> _bamStreams;
};
