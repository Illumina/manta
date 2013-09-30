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

#include "options/ReadScannerOptionsParser.hh"


boost::program_options::options_description
getOptionsDescription(ReadScannerOptions& opt)
{
    namespace po = boost::program_options;
    po::options_description desc("read-scanner");
    desc.add_options()
    ("min-candidate-sv-size", po::value(&opt.minCandidateVariantSize)->default_value(opt.minCandidateVariantSize),
     "Indels below this size will not be discovered or reported as candidates")
    ("min-mapq", po::value(&opt.minMapq)->default_value(opt.minMapq),
     "Reads with MAPQ less than this value will be ignored")
    ("edge-prob", po::value(&opt.breakendEdgeTrimProb)->default_value(opt.breakendEdgeTrimProb),
     "Breakend range associated with each read will have this probability trimmed from each edge")
    ;

    return desc;
}



bool
parseOptions(
    const boost::program_options::variables_map& /*vm*/,
    ReadScannerOptions& opt,
    std::string& errorMsg)
{
    errorMsg.clear();
    if ((opt.breakendEdgeTrimProb <= 0) || (opt.breakendEdgeTrimProb >= 1.0))
    {
        errorMsg="edge-prob argument is restricted to (0,1)";
    }

    return (! errorMsg.empty());

}
