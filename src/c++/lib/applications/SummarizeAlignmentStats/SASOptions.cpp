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

#include "SASOptions.hpp"

#include "blt_util/log.hpp"
#include "options/optionsUtil.hpp"

#include "boost/program_options.hpp"

#include <iostream>
#include "../../common/ProgramUtil.hpp"

static void usage(
    std::ostream&                                      os,
    const illumina::Program&                           prog,
    const boost::program_options::options_description& visible,
    const char*                                        msg = nullptr)
{
  usage(os, prog, visible, "summarize fragment size stats", " [ > output]", msg);
}

static void checkStandardizeUsageFile(
    std::ostream&                                      os,
    const illumina::Program&                           prog,
    const boost::program_options::options_description& visible,
    std::string&                                       filename,
    const char*                                        fileLabel)
{
  std::string errorMsg;
  if (checkAndStandardizeRequiredInputFilePath(filename, fileLabel, errorMsg)) {
    usage(os, prog, visible, errorMsg.c_str());
  }
}

void parseSASOptions(const illumina::Program& prog, int argc, char* argv[], SASOptions& opt)
{
  namespace po = boost::program_options;
  po::options_description req("configuration");
  // clang-format off
  req.add_options()
  ("align-stats", po::value(&opt.statsFilename),
   "pre-computed alignment statistics for the input alignment files (required)")
  ("output-file", po::value(&opt.outputFilename),
   "write fragment summary stats to filename (default: stdout)");
  ;
  // clang-format on

  po::options_description help("help");
  help.add_options()("help,h", "print this message");

  po::options_description visible("options");
  visible.add(req).add(help);

  bool              po_parse_fail(false);
  po::variables_map vm;
  try {
    po::store(
        po::parse_command_line(
            argc, argv, visible, po::command_line_style::unix_style ^ po::command_line_style::allow_short),
        vm);
    po::notify(vm);
  } catch (const boost::program_options::error&
               e)  // todo:: find out what is the more specific exception class thrown by program options
  {
    log_os << "\nERROR: Exception thrown by option parser: " << e.what() << "\n";
    po_parse_fail = true;
  }

  if ((argc <= 1) || (vm.count("help")) || po_parse_fail) {
    usage(log_os, prog, visible);
  }

  checkStandardizeUsageFile(log_os, prog, visible, opt.statsFilename, "alignment statistics");
}
