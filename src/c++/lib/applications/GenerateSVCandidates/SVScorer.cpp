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

#include "SVScorer.hh"

#include "blt_util/bam_streamer.hh"
#include "blt_util/log.hh"
#include "common/Exceptions.hh"
#include "manta/ReadGroupStatsSet.hh"

#include "boost/foreach.hpp"

#include <iostream>



SVScorer::
SVScorer(const GSCOptions& opt) :
    _somaticOpt(opt.somaticOpt),
    _readScanner(opt.scanOpt,opt.statsFilename,opt.alignmentFilename)
{
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
SVScorer::
scoreSomaticSV(
    const SVCandidateData& svData,
    const unsigned svIndex,
    const SVCandidate& sv,
    SomaticSVScoreInfo& ssInfo)
{

}




