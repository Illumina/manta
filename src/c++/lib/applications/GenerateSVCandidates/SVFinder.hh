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
#include "EdgeInfo.hh"

#include "blt_util/bam_streamer.hh"
#include "manta/SVCandidate.hh"
#include "manta/SVLocusSet.hh"

#include "boost/shared_ptr.hpp"

#include <vector>


struct SVFinder {

    SVFinder(const GSCOptions& opt);

    const SVLocusSet&
    getSet() const
    {
        return _set;
    }

    void
    findSVCandidates(
            const EdgeInfo& edge,
            std::vector<SVCandidate>& svs);

private:
    SVLocusSet _set;

    typedef boost::shared_ptr<bam_streamer> streamPtr;
    std::vector<streamPtr> _bamStreams;
};
