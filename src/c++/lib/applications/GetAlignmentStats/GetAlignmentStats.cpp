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

#include "GetAlignmentStats.hh"
#include "AlignmentStatsOptions.hh"

#include "blt_util/log.hh"
#include "manta/ReadGroupStatsSet.hh"

#include "boost/foreach.hpp"

#include <cstdlib>

#include <iostream>



static
void
runAlignmentStats(const AlignmentStatsOptions& opt) {

    // calculate mean, median and standard deviation of the insert
    // size for each bam file
    ReadGroupStatsSet rstats;

    if (opt.alignmentFilename.empty()) {
        log_os << "ERROR: No input files specified.\n";
        exit(EXIT_FAILURE);
    }

    BOOST_FOREACH(const std::string& file, opt.alignmentFilename) {
        rstats.setStats(file,ReadGroupStats(file));
    }

    std::ostream& statfp(std::cout);
    rstats.store(statfp);
}



void
GetAlignmentStats::
runInternal(int argc, char* argv[]) const {

    AlignmentStatsOptions opt;

    parseAlignmentStatsOptions(*this,argc,argv,opt);
    runAlignmentStats(opt);
}
