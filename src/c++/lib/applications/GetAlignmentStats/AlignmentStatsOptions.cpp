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

/// \file

#include "AlignmentStatsOptions.hh"

#include "blt_util/log.hh"

#include "boost/program_options.hpp"

#include <iostream>



void
parseAlignmentStatsOptions(const manta::Program& prog,
                           int argc, char* argv[],
                           AlignmentStatsOptions& opt) {

    namespace po = boost::program_options;
    po::options_description req("configuration");
    req.add_options()
    ("align-file", po::value<std::vector<std::string> >(&opt.alignmentFilename),
     "alignment file in bam format (may be specified multiple times)")
    ("output-file", po::value<std::string>(&opt.outputFilename),
     "write stats to filename (default: stdout)");

    po::options_description help("help");
    help.add_options()
    ("help,h","print this message");

    po::options_description visible("options");
    visible.add(req).add(help);

    bool po_parse_fail(false);
    po::variables_map vm;
    try {
        po::store(po::parse_command_line(argc, argv, visible,
                                         po::command_line_style::unix_style ^ po::command_line_style::allow_short), vm);
        po::notify(vm);
    } catch (const boost::program_options::error& e) { // todo:: find out what is the more specific exception class thrown by program options
        log_os << "\nERROR: Exception thrown by option parser: " << e.what() << "\n";
        po_parse_fail=true;
    }

    if ((argc<=1) || (vm.count("help")) || po_parse_fail) {
        log_os << "\n" << prog.name() << ": get statistics for SV-calling from alignment files\n\n";
        log_os << "version: " << prog.version() << "\n\n";
        log_os << "usage: " << prog.name() << " [options] > stats\n\n";
        log_os << visible << "\n";
        exit(EXIT_FAILURE);
    }
}

