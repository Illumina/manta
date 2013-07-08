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

#pragma once

#include "GSCOptions.hh"
#include "manta/EdgeInfo.hh"

#include "blt_util/bam_streamer.hh"
#include "manta/SVCandidate.hh"
#include "manta/SVCandidateData.hh"
#include "manta/SVLocusScanner.hh"
#include "manta/SVLocusSet.hh"

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
        const EdgeInfo& edge,
        SVCandidateData& svData,
        std::vector<SVCandidate>& svs);

    void
    checkResult(
        const SVCandidateData& svData,
        const std::vector<SVCandidate>& svs) const;

private:

    void
    addSVNodeData(
        const SVLocus& locus,
        const NodeIndexType node1,
        const NodeIndexType node2,
        SVCandidateData& svData);


    void
    getCandidatesFromData(
        SVCandidateData& svData,
        std::vector<SVCandidate>& svs);

    const ReadScannerOptions _scanOpt;
    SVLocusSet _set;
    SVLocusScanner _readScanner;

    typedef boost::shared_ptr<bam_streamer> streamPtr;
    std::vector<streamPtr> _bamStreams;
};
