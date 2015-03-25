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

#pragma once

#include "options/AlignmentFileOptions.hh"

#include "boost/program_options.hpp"

#include <string>


boost::program_options::options_description
getOptionsDescription(
    AlignmentFileOptions& opt);


bool
parseOptions(
    const boost::program_options::variables_map& vm,
    AlignmentFileOptions& opt,
    std::string& errorMsg);
