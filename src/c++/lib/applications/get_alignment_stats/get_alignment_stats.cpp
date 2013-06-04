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

#include "get_alignment_stats.hh"
#include "alignment_stats_options.hh"

#include "blt_util/log.hh"
#include "manta_common/read_group_stats_set.hh"

#include "boost/foreach.hpp"

#include <cstdlib>

#include <iostream>



static
void
run_alignment_stats(const alignment_stats_options& opt) {

    // calculate mean, median and standard deviation of the insert
    // size for each bam file
    read_group_stats_set rstats;

    if (opt.alignment_filename.empty()) {
        log_os << "ERROR: No input files specified.\n";
        exit(EXIT_FAILURE);
    }

    BOOST_FOREACH(const std::string& file, opt.alignment_filename) {
        rstats.set_stats(file,read_group_stats(file));
    }

    std::ostream& statfp(std::cout);
    rstats.store(statfp);
}



void
get_alignment_stats::
run_internal(int argc, char* argv[]) const {

    alignment_stats_options opt;

    parse_alignment_stats_options(*this,argc,argv,opt);
    run_alignment_stats(opt);
}
