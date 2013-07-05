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

#pragma once

#include "boost/program_options.hpp"


struct ReadScannerOptions
{

    ReadScannerOptions() :
        minMapq(15),
        breakendEdgeTrimProb(0.25)
    {}

    unsigned minMapq;

    /// report breakend regions with x prob regions removed from each edge
    float breakendEdgeTrimProb;
};


boost::program_options::options_description
getOptionsDescription(ReadScannerOptions& opt);

