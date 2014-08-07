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


struct ReadScannerOptions
{
    ReadScannerOptions() {}

    /// standard MAPQ filter applied during locus generation and some/not all subsequent steps
    unsigned minMapq = 15;

    /// a second, lower mapq threshold used only during somatic calling to disprove a
    /// somatic candidate using weak normal sample evidence
    unsigned minTier2Mapq = 5;

    /// report breakend regions with x prob regions removed from each edge
    float breakendEdgeTrimProb = 0.25;

    /// report breakend regions with x prob regions removed from each edge
    /// used only for 'large-scale' events.
    float largeScaleEventBreakendEdgeTrimProb = 0.1;

    /// report a pair as "proper pair" if fragment size is within x prob region removed from each edge
    float properPairTrimProb = 0.01;

    /// add a pair to the evidence pool if frag size is within x prob region removed from each edge
    float evidenceTrimProb = 0.15;

    /// fragment length to search upstream of a breakend for shadow read support
    float shadowSearchRangeProb = 0.05;

    /// multiplier for fragment length
    float shadowSearchRangeFactor = 1.2;

    ///< ignore indels smaller than this when building graph:
    unsigned minCandidateVariantSize = 10;

    /// whenever a breakend is predicted from a read pair junction, the predicted breakend
    /// range should be no smaller than this:
    unsigned minPairBreakendSize = 40;

    /// whenever a breakend is predicted from an individual read split (ie. non-assembled),
    /// set the predicted breakend size to this fraction of the
    /// event size (modified by the min and max limits below)
    float splitBreakendSizeFraction = 0.1;

    /// whenever a breakend is predicted from an individual read split (ie. non-assembled),
    /// the predicted breakend range should be no larger than this:
    unsigned maxSplitBreakendSize = 100;

    /// whenever a breakend is predicted from an individual read split (ie. non-assembled),
    /// the predicted breakend range should be no smaller than this:
    unsigned minSplitBreakendSize = 10;

    /// Semi-aligned regions (including soft-clipped) need to be at least this long to be included as SV evidence
    ///
    unsigned minSemiAlignedMismatchLen = 8;

    /// Minimal length of a deletion / insertion SV candidate in RNA data
    unsigned minRNALength = 100000;

    /// Accept semi-aligned reads with at least this hypothesis score
    /// different for graph and candidate generation
    double minSemiAlignedScoreGraph = 180.0;
    double minSemiAlignedScoreCandidates = 180.0;

    ///< min MAPQ for shadow mate used to build SV adjacency graph
    unsigned minSingletonMapqGraph = 30;

    /// min MAPQ for shadow mate used for candidate assembly and scoring
    unsigned minSingletonMapqCandidates = 15;

    /// typically set true for RNA-Seq analysis, where proper-pair is used to signal intron-spanning pairs
    bool isIgnoreAnomProperPair = false;

    /// the maximum depth at which input reads are considered in graph creation/assembly, etc.
    /// (when avg chrom depths are provided)
    float maxDepthFactor = 12;

    /// the maximum depth for a whole locus for remote read retrieval (ie. MAPQ0 chimera mates retrieved for large insertion assembly)
    float maxDepthFactorRemoteReads = 7;
};

