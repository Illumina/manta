//
// Manta - Structural Variant and Indel Caller
// Copyright (c) 2013-2025 Illumina, Inc.
//
// This program is licensed under the terms of the Polyform strict license
//
// ***As far as the law allows, the software comes as is, without
// any warranty or condition, and the licensor will not be liable
// to you for any damages arising out of these terms or the use
// or nature of the software, under any kind of legal claim.***
//
// You should have received a copy of the PolyForm Strict License 1.0.0
// along with this program.  If not, see <https://polyformproject.org/licenses/strict/1.0.0>.
//
//

/// \file
/// \author Chris Saunders
///

#include "ProgramUtil.hpp"

#include <iostream>

void usage(
    std::ostream&                                      os,
    const illumina::Program&                           prog,
    const boost::program_options::options_description& visible,
    const char*                                        desc,
    const char*                                        afteropts,
    const char*                                        msg)
{
  os << "\n" << prog.name() << ": " << desc << "\n\n";
  os << "version: " << prog.version() << "\n";
  os << "compiler: " << prog.compiler() << "\n";
  os << "build-time: " << prog.buildTime() << "\n\n";
  os << "usage: " << prog.name() << " [options]" << afteropts << "\n\n";
  os << visible << "\n\n";

  if (nullptr != msg) {
    os << msg << "\n\n";
  }
  exit(2);
}
