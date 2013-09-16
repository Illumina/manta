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

#pragma once

#include "EdgeOptions.hh"

#include "boost/program_options.hpp"

#include <string>


boost::program_options::options_description
getOptionsDescription(
    EdgeOptions& opt);


/// additional parsing beyond default boost::program_options behaviors,
/// and arg validation
///
/// \return is valid
bool
parseOptions(
    const boost::program_options::variables_map& vm,
    EdgeOptions& opt,
    std::string& errorMsg);
