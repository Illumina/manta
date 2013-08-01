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

#include "GSCOptions.hh"

#include "blt_util/log.hh"

#include "boost/filesystem.hpp"
#include "boost/foreach.hpp"
#include "boost/program_options.hpp"

#include <iostream>
#include <set>



static
void
usage(
    std::ostream& os,
    const manta::Program& prog,
    const boost::program_options::options_description& visible,
    const char* msg = NULL)
{
    os << "\n" << prog.name() << ": call candidates from an SV Locus graph\n\n";
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
parseGSCOptions(const manta::Program& prog,
                int argc, char* argv[],
                GSCOptions& opt)
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
    ("graph-file", po::value(&opt.graphFilename),
     "sv locus graph file (required)")
    ("align-stats", po::value(&opt.statsFilename),
     "pre-computed alignment statistics for the input alignment files (required)")
    ("chrom-depth", po::value(&opt.chromDepthFilename),
     "average depth estimate for each chromosome")
    ("ref", po::value(&opt.referenceFilename),
     "fasta reference sequence (required)")
    ("candidate-output-file", po::value(&opt.candidateOutputFilename),
     "Write SV candidates to file (required)")
    ("somatic-output-file", po::value(&opt.somaticOutputFilename),
     "Write somatic SV candidates to file (at least one tumor and non-tumor alignment file must be specified)")
    ("bin-count", po::value(&opt.binCount)->default_value(opt.binCount),
     "Specify how many bins the SV candidate problem should be divided into, where bin-index can be used to specify which bin to solve")
    ("bin-index", po::value(&opt.binIndex)->default_value(opt.binIndex),
     "specify which bin to solve when the SV candidate problem is subdivided into bins. Value must bin in [0,bin-count)")
    ("locus-index", po::value(&opt.locusIndex),
     "Instead of solving for all SV candidates in a bin, solve for candidates of a particular locus. If specified then bin-index is ignored.")
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

    opt.isLocusIndex=vm.count("locus-index");

    {
        // paste together tumor and normal:
        opt.alignmentFilename = normalAlignmentFilename;
        opt.alignmentFilename.insert(opt.alignmentFilename.end(),
                                     tumorAlignmentFilename.begin(),
                                     tumorAlignmentFilename.end());
        opt.isAlignmentTumor.clear();
        opt.isAlignmentTumor.resize(normalAlignmentFilename.size(),false);
        opt.isAlignmentTumor.resize(opt.alignmentFilename.size(),true);
    }

    // fast check of config state:
    if (opt.binCount < 1)
    {
        usage(log_os,prog,visible,"bin-count must be 1 or greater");
    }
    if (opt.binIndex >= opt.binCount)
    {
        usage(log_os,prog,visible,"bin-index must be in range [0,bin-count)");
    }
    if (opt.alignmentFilename.empty())
    {
        usage(log_os,prog,visible,"Must specify at least one input alignment file");
    }
    if (opt.statsFilename.empty())
    {
        usage(log_os,prog,visible,"Need the alignment stats file");
    }
    if (opt.referenceFilename.empty())
    {
    	usage(log_os,prog,visible,"Need the FASTA reference file");
    }
    {
        // check that alignment files exist, and names do not repeat
        std::set<std::string> nameCheck;
        BOOST_FOREACH(std::string& afile, opt.alignmentFilename)
        {
            checkStandardizeUsageFile(log_os,prog,visible,afile,"alignment file");
            if (nameCheck.count(afile))
            {
                std::ostringstream oss;
                oss << "Repeated alignment filename: " << afile << "\n";
                usage(log_os,prog,visible,oss.str().c_str());
            }
            nameCheck.insert(afile);
        }
    }
    checkStandardizeUsageFile(log_os,prog,visible,opt.graphFilename,"SV locus graph");
    checkStandardizeUsageFile(log_os,prog,visible,opt.referenceFilename,"reference fasta");
    checkStandardizeUsageFile(log_os,prog,visible,opt.statsFilename,"alignment statistics");

    if (! opt.chromDepthFilename.empty())
    {
        checkStandardizeUsageFile(log_os,prog,visible,opt.chromDepthFilename,"chromosome depth");
    }
    if (opt.candidateOutputFilename.empty())
    {
        usage(log_os,prog,visible,"Must specify candidate output file");
    }

    if (! opt.somaticOutputFilename.empty())
    {
        if (normalAlignmentFilename.empty() || tumorAlignmentFilename.empty())
        {
            usage(log_os,prog,visible,"Must specify at least one tumor and non-tumor alignment file for somatic output");
        }
    }
}

