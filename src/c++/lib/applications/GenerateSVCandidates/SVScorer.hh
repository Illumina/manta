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

#include "blt_util/bam_streamer.hh"
#include "manta/SVCandidate.hh"
#include "manta/SVCandidateData.hh"
#include "manta/SVLocusScanner.hh"
#include "manta/SomaticSVScoreInfo.hh"

#include "boost/shared_ptr.hpp"

#include <vector>


struct SVScorer
{
    SVScorer(const GSCOptions& opt);

    void
    scoreSomaticSV(
        const SVCandidateData& svData,
        const unsigned svIndex,
        const SVCandidate& sv,
        SomaticSVScoreInfo& ssInfo);

private:
    const SomaticCallOptions _somaticOpt;
    SVLocusScanner _readScanner;

    typedef boost::shared_ptr<bam_streamer> streamPtr;
    std::vector<streamPtr> _bamStreams;
};
