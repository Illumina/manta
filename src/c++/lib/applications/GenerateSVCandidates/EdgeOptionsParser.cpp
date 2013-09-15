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
// <https://github.com/sequencing/licenses/>
//

///
/// \author Chris Saunders
///

#include "EdgeOptionsParser.hh"



boost::program_options::options_description
getOptionsDescription(
    EdgeOptions& opt)
{
    namespace po = boost::program_options;

    po::options_description optdesc("edge-selection");
    optdesc.add_options()
    ("bin-count", po::value(&opt.binCount)->default_value(opt.binCount),
     "Specify how many bins the SV candidate problem should be divided into, where bin-index can be used to specify which bin to solve")
    ("bin-index", po::value(&opt.binIndex)->default_value(opt.binIndex),
     "specify which bin to solve when the SV candidate problem is subdivided into bins. Value must bin in [0,bin-count)")
    ("locus-index", po::value(&opt.locusIndex),
     "Instead of solving for all SV candidates in a bin, solve for candidates of a particular locus. If specified then bin-index is ignored.")
    ;
    return optdesc;
}



bool
parseOptions(
    const boost::program_options::variables_map& vm,
    EdgeOptions& opt,
    std::string& errorMsg)
{
    errorMsg.clear();

    opt.isLocusIndex=vm.count("locus-index");

    // fast check of config state:
    if (opt.binCount < 1)
    {
        errorMsg="bin-count must be 1 or greater";
    }
    else if (opt.binIndex >= opt.binCount)
    {
        errorMsg="bin-index must be in range [0,bin-count)";
    }

    return (! errorMsg.empty());
}
