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

/// \author Chris Saunders
///

#include "ProgramUtil.hh"

#include <iostream>



void
usage(
    std::ostream& os,
    const manta::Program& prog,
    const boost::program_options::options_description& visible,
    const char* desc,
    const char* afteropts,
    const char* msg)
{
    os << "\n" << prog.name() << ": " << desc << "\n\n";
    os << "version: " << prog.version() << "\n\n";
    os << "compiler: " << prog.compiler() << "\n\n";
    os << "build-time: " << prog.buildTime() << "\n\n";
    os << "usage: " << prog.name() << " [options]" << afteropts << "\n\n";
    os << visible << "\n\n";

    if (nullptr != msg)
    {
        os << msg << "\n\n";
    }
    exit(2);
}
