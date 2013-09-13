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

#include "boost/program_options.hpp"


struct ReadScannerOptions
{
    ReadScannerOptions() :
        minMapq(15),
        breakendEdgeTrimProb(0.25),
        properPairTrimProb(0.01),
        evidenceTrimProb(0.15),
        minCandidateIndelSize(10),
        minPairBreakendSize(40),
        splitBreakendSizeFraction(0.1),
        maxSplitBreakendSize(100),
        minSplitBreakendSize(10)
    {}

    unsigned minMapq;

    /// report breakend regions with x prob regions removed from each edge
    float breakendEdgeTrimProb;

    /// report a pair as "proper pair" if fragment size is within x prob region removed from each edge
    float properPairTrimProb;

    /// add a pair to the evidence pool if frag size is within x prob region removed from each edge
    float evidenceTrimProb;

    /// ignore indels smaller than this when building graph:
    unsigned minCandidateIndelSize;

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
};


boost::program_options::options_description
getOptionsDescription(ReadScannerOptions& opt);
