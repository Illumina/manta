// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Manta
// Copyright (c) 2013-2015 Illumina, Inc.
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

#include "SSLOptions.hh"

#include "blt_util/log.hh"

#include "boost/filesystem.hpp"
#include "boost/program_options.hpp"

#include <iostream>



static
void
usage(
    std::ostream& os,
    const manta::Program& prog,
    const boost::program_options::options_description& visible,
    const char* msg = nullptr)
{
    os << "\n" << prog.name() << ": write graph summary stats to stdout (tsv format)\n\n";
    os << "version: " << prog.version() << "\n\n";
    os << "usage: " << prog.name() << " [options] > graph_summary\n\n";
    os << visible << "\n\n";

    if (NULL != msg)
    {
        os << msg << "\n\n";
    }
    exit(2);
}



void
parseSSLOptions(
    const manta::Program& prog,
    int argc, char* argv[],
    SSLOptions& opt)
{
    namespace po = boost::program_options;
    po::options_description req("configuration");
    req.add_options()
    ("graph-file", po::value(&opt.graphFilename),
     "sv locus graph file")
    ("global", po::value(&opt.isGlobalStats)->zero_tokens(),
     "provide global stats on full graph (default output is per-locus stats)")
    ;

    po::options_description help("help");
    help.add_options()
    ("help,h","print this message");

    po::options_description visible("options");
    visible.add(req).add(help);

    bool po_parse_fail(false);
    po::variables_map vm;
    try
    {
        po::store(po::parse_command_line(argc, argv, visible,
                                         po::command_line_style::unix_style ^ po::command_line_style::allow_short), vm);
        po::notify(vm);
    }
    catch (const boost::program_options::error& e)     // todo:: find out what is the more specific exception class thrown by program options
    {
        log_os << "\nERROR: Exception thrown by option parser: " << e.what() << "\n";
        po_parse_fail=true;
    }

    if ((argc<=1) || (vm.count("help")) || po_parse_fail)
    {
        usage(log_os,prog,visible);
    }

    // fast check of config state:
    if (opt.graphFilename.empty())
    {
        usage(log_os,prog,visible,"Must specify sv locus graph file");
    }
    if (! boost::filesystem::exists(opt.graphFilename))
    {
        usage(log_os,prog,visible,"SV locus graph file does not exist");
    }
}

