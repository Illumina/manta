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

#include "options/SVLocusSetOptionsParser.hh"


boost::program_options::options_description
getOptionsDescription(SVLocusSetOptions& opt)
{
    namespace po = boost::program_options;
    po::options_description desc("sv-locus-graph");
    desc.add_options()
    ("min-edge-observations", po::value(&opt.minMergeEdgeObservations)->default_value(opt.minMergeEdgeObservations),
     "Minimum number of supporting observations required to retain a graph edge")
    ;

    return desc;
}



bool
parseOptions(
    const boost::program_options::variables_map& /*vm*/,
    SVLocusSetOptions& /*opt*/,
    std::string& errorMsg)
{
    errorMsg.clear();
#if 0
    if ((opt.breakendEdgeTrimProb <= 0) || (opt.breakendEdgeTrimProb >= 1.0))
    {
        errorMsg="edge-prob argument is restricted to (0,1)";
    }
#endif
    return (! errorMsg.empty());

}
