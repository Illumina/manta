// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Manta
// Copyright (c) 2013-2014 Illumina, Inc.
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

#include "ESLOptions.hh"

#include "blt_util/log.hh"
#include "options/AlignmentFileOptionsParser.hh"
#include "options/ReadScannerOptionsParser.hh"
#include "options/SVLocusSetOptionsParser.hh"
#include "options/optionsUtil.hh"

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
    os << "\n" << prog.name() << ": partition sv evidence regions\n\n";
    os << "version: " << prog.version() << "\n\n";
    os << "usage: " << prog.name() << " [options]\n\n";
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
    std::string errorMsg;
    if ( checkStandardizeInputFile(filename, fileLabel, errorMsg))
    {
        usage(os,prog,visible,errorMsg.c_str());
    }
}



void
parseESLOptions(
    const manta::Program& prog,
    int argc, char* argv[],
    ESLOptions& opt)
{
    namespace po = boost::program_options;
    po::options_description req("configuration");
    req.add_options()
    ("output-file", po::value(&opt.outputFilename),
     "write SV Locus graph to file (required)")
    ("ref", po::value(&opt.referenceFilename),
     "fasta reference sequence (required)")
    ("align-stats", po::value(&opt.statsFilename),
     "pre-computed alignment statistics for the input alignment files (required)")
    ("chrom-depth", po::value(&opt.chromDepthFilename),
     "average depth estimate for each chromosome")
    ("truth-vcf", po::value(&opt.truthVcfFilename),
     "optional truth VCF file (for testing)")
    ("region", po::value(&opt.region),
     "samtools formatted region, eg. 'chr1:20-30'");

    po::options_description alignDesc(getOptionsDescription(opt.alignFileOpt));
    po::options_description scanDesc(getOptionsDescription(opt.scanOpt));
    po::options_description graphDesc(getOptionsDescription(opt.graphOpt));

    po::options_description help("help");
    help.add_options()
    ("help,h","print this message");

    po::options_description visible("options");
    visible.add(alignDesc).add(scanDesc).add(graphDesc).add(req).add(help);

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

    std::string errorMsg;
    if (parseOptions(vm, opt.alignFileOpt, errorMsg))
    {
        usage(log_os,prog,visible,errorMsg.c_str());
    }
    else if (parseOptions(vm, opt.scanOpt, errorMsg))
    {
        usage(log_os,prog,visible,errorMsg.c_str());
    }
    else if (opt.outputFilename.empty())
    {
        usage(log_os,prog,visible,"Must specify a graph output file");
    }
    else if (opt.referenceFilename.empty())
    {
        usage(log_os,prog,visible,"Need the FASTA reference file");
    }
    else if (opt.region.empty())
    {
        usage(log_os,prog,visible,"Need the samtools formatted region");
    }

    checkStandardizeUsageFile(log_os,prog,visible,opt.statsFilename,"alignment statistics");
    checkStandardizeUsageFile(log_os,prog,visible,opt.referenceFilename,"reference fasta");

    if (! opt.chromDepthFilename.empty())
    {
        checkStandardizeUsageFile(log_os,prog,visible,opt.chromDepthFilename,"chromosome depth");
    }
}
