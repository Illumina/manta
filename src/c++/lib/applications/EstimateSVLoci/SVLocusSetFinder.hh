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

#pragma once

#include "blt_util/bam_record.hh"
#include "manta/ReadGroupStatsSet.hh"
#include "manta/SVLocus.hh"

#include <string>
#include <vector>


// estimate an SVLocusSet
//
struct SVLocusSetFinder {

    SVLocusSetFinder(
            const std::string& statsFilename,
            const std::string& region)
    {
        // pull in insert stats:
        _rss.load(statsFilename.c_str());
    }

    void
    update(const bam_record& read,
           const unsigned defaultReadGroupIndex) {}

private:
    ReadGroupStatsSet _rss;
    std::vector<SVLocus> _loci;
};

