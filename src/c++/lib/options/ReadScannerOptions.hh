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


struct ReadScannerOptions
{
    ReadScannerOptions() :
        minMapq(15),
        breakendEdgeTrimProb(0.25),
        properPairTrimProb(0.01),
        evidenceTrimProb(0.15),
        minCandidateVariantSize(10),
        minPairBreakendSize(40),
        splitBreakendSizeFraction(0.1),
        maxSplitBreakendSize(100),
        minSplitBreakendSize(10),
        minSoftClipLen(8),
        minSemiAlignedScoreGraph(50.0),
        minSemiAlignedScoreCandidates(30.0),
        minSingletonMapqGraph(40),
        minSingletonMapqCandidates(20),
        isIgnoreAnomProperPair(false)
    {}

    unsigned minMapq;

    /// report breakend regions with x prob regions removed from each edge
    float breakendEdgeTrimProb;

    /// report a pair as "proper pair" if fragment size is within x prob region removed from each edge
    float properPairTrimProb;

    /// add a pair to the evidence pool if frag size is within x prob region removed from each edge
    float evidenceTrimProb;

    /// ignore indels smaller than this when building graph:
    unsigned minCandidateVariantSize;

    // whenever a breakend is predicted from a read pair junction, the predicted breakend range should be no
    // smaller than this:
    unsigned minPairBreakendSize;

    // whenever a breakend is predicted from an individual read split (ie. non-assembled), set the predicted breakend size to this fraction of the
    // event size (modified by the min and max limits below)
    float splitBreakendSizeFraction;

    // whenever a breakend is predicted from an individual read split (ie. non-assembled), the predicted breakend range should be no
    // larger than this:
    unsigned maxSplitBreakendSize;

    // whenever a breakend is predicted from an individual read split (ie. non-assembled), the predicted breakend range should be no
    // smaller than this:
    unsigned minSplitBreakendSize;

    // Soft clipped read ends must be of at least this length to be entered as small SV evidence
    unsigned minSoftClipLen;

    // Accept semi-aligned reads with at least this hypothesis score, different for graph and candidate generation
    double minSemiAlignedScoreGraph;
    double minSemiAlignedScoreCandidates;

    // We want only shadows with a good singleton mapq, but use again different thresholds for graph and candidate generation
    unsigned minSingletonMapqGraph;
    unsigned minSingletonMapqCandidates;

    bool isIgnoreAnomProperPair; ///< typically set true for RNA-Seq analysis, where proper-pair is used to signal intron-spanning pairs
};
