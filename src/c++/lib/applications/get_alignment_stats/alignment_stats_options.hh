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

#include "manta_common/program.hh"

#include <string>
#include <vector>



struct alignment_stats_options {

    std::vector<std::string> alignment_filename;
};


void
parse_alignment_stats_options(const manta::program& prog,
                              int argc, char* argv[],
                              alignment_stats_options& opt);
