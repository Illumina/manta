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

#include "boost/filesystem.hpp"
#include "boost/foreach.hpp"
#include "boost/program_options.hpp"

#include <iostream>



static
void
usage(
    std::ostream& os,
    const manta::Program& prog,
    const boost::program_options::options_description& visible,
    const char* msg = NULL)
{
    os << "\n" << prog.name() << ": get statistics for SV-calling from alignment files\n\n";
    os << "version: " << prog.version() << "\n\n";
    os << "usage: " << prog.name() << " [options] > stats\n\n";
    os << visible << "\n\n";

    if (NULL != msg)
    {
        os << msg << "\n\n";
    }
    exit(2);
}



static
void
checkStandardizeUsageFile(
    std::ostream& os,
    const manta::Program& prog,
    const boost::program_options::options_description& visible,
    std::string& filename,
    const char* fileLabel)
{
    if (filename.empty())
    {
        std::ostringstream oss;
        oss << "Must specify " << fileLabel << " file";
        usage(os,prog,visible,oss.str().c_str());
    }
    if (! boost::filesystem::exists(filename))
    {
        std::ostringstream oss;
        oss << "Can't find " << fileLabel << " file '" << filename << "'";
        usage(os,prog,visible,oss.str().c_str());
    }
    filename = boost::filesystem::canonical(filename).string();
}



void
parseAlignmentStatsOptions(const manta::Program& prog,
                           int argc, char* argv[],
                           AlignmentStatsOptions& opt)
{
    std::vector<std::string> normalAlignmentFilename;
    std::vector<std::string> tumorAlignmentFilename;

    namespace po = boost::program_options;
    po::options_description req("configuration");
    req.add_options()
    ("align-file", po::value(&normalAlignmentFilename),
     "alignment file in bam format (may be specified multiple times, assumed to be non-tumor if tumor file(s) provided)")
    ("tumor-align-file", po::value(&tumorAlignmentFilename),
     "tumor sample alignment file in bam format (may be specified multiple times)")
    ("output-file", po::value(&opt.outputFilename),
     "write stats to filename (default: stdout)");

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
        log_os << "\n" << prog.name() << ": get statistics for SV-calling from alignment files\n\n";
        log_os << "version: " << prog.version() << "\n\n";
        log_os << "usage: " << prog.name() << " [options] > stats\n\n";
        log_os << visible << "\n";
        exit(EXIT_FAILURE);
    }

    {
        // paste together tumor and normal:
        opt.alignmentFilename = normalAlignmentFilename;
        opt.alignmentFilename.insert(opt.alignmentFilename.end(),
                                     tumorAlignmentFilename.begin(),
                                     tumorAlignmentFilename.end());
    }


    // fast check of config state:
    if (opt.alignmentFilename.empty())
    {
        usage(log_os,prog,visible,"Must specify at least one input alignment file");
    }
    {
        // check that alignment files exist, and names do not repeat
        std::set<std::string> nameCheck;
        BOOST_FOREACH(std::string& afile, opt.alignmentFilename)
        {
            checkStandardizeUsageFile(log_os,prog,visible,afile,"alignment file");
            if(nameCheck.count(afile))
            {
                std::ostringstream oss;
                oss << "Repeated alignment filename: " << afile << "\n";
                usage(log_os,prog,visible,oss.str().c_str());
            }
            nameCheck.insert(afile);
        }
    }
}

