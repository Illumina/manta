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

#pragma once

#include "GSCOptions.hh"

#include "blt_util/bam_streamer.hh"
#include "manta/ChromDepthFilterUtil.hh"
#include "manta/SVCandidate.hh"
#include "manta/SVCandidateSetData.hh"
#include "manta/SVLocusScanner.hh"
#include "svgraph/EdgeInfo.hh"
#include "svgraph/SVLocusSet.hh"

#include <vector>


struct SVFinder
{
    SVFinder(const GSCOptions& opt);

    ~SVFinder();

    const SVLocusSet&
    getSet() const
    {
        return _set;
    }

    void
    findCandidateSV(
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
        const bam_header_info& bamHeader,
        const SVLocus& locus,
        const NodeIndexType node1,
        const NodeIndexType node2,
        const GenomeInterval& searchInterval,
        const reference_contig_segment& refSeq,
        const bool isNode1,
        const bool isSomatic,
        SVCandidateSetData& svData,
        TruthTracker& truthTracker);

    void
    assignPairObservationsToSVCandidates(
        const SVLocusNode& node1,
        const SVLocusNode& node2,
        const std::vector<SVObservation>& readCandidates,
        const bool isExpandSVCandidateSet,
        SVCandidateSetReadPair& pair,
        std::vector<SVCandidate>& svs);

    /// we either process the pair to discover new SVs and expand existing SVs,
    /// or we go through and add pairs to existing SVs without expansion
    ///
    void
    processReadPair(
        const SVLocusNode& node1,
        const SVLocusNode& node2,
        const bam_header_info& bamHeader,
        const reference_contig_segment& refSeq1,
        const reference_contig_segment& refSeq2,
        const unsigned bamIndex,
        const bool isExpandSVCandidateSet,
        std::vector<SVCandidate>& svs,
        TruthTracker& truthTracker,
        SVCandidateSetReadPair& pair);

    void
    getCandidatesFromData(
        const SVLocusNode& node1,
        const SVLocusNode& node2,
        const bam_header_info& bamHeader,
        const reference_contig_segment& refSeq1,
        const reference_contig_segment& refSeq2,
        const bool isSomatic,
        SVCandidateSetData& svData,
        std::vector<SVCandidate>& svs,
        TruthTracker& truthTracker);

    const ChromDepthFilterUtil&
    dFilter() const
    {
        return *(_dFilterPtr);
    }

    const ReadScannerOptions _scanOpt;
    const std::vector<bool> _isAlignmentTumor;
    SVLocusSet _set;
    std::unique_ptr<ChromDepthFilterUtil> _dFilterPtr;
    SVLocusScanner _readScanner;

    const std::string _referenceFilename;

    const bool _isRNA;

    typedef std::shared_ptr<bam_streamer> streamPtr;
    std::vector<streamPtr> _bamStreams;

    /// this is only here as syscall cache:
    std::vector<SVObservation> _readCandidates;

    /// throwaway stats tracker...
    SampleEvidenceCounts _eCounts;
};
