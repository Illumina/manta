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

#include "MergeAlignmentStatsOptions.hpp"
#include "blt_util/log.hpp"
#include "options/AlignmentFileOptionsParser.hpp"
#include "options/optionsUtil.hpp"

#include "boost/program_options.hpp"

#include <iostream>
#include "../../common/ProgramUtil.hpp"

typedef std::vector<std::string> files_t;

static void usage(
    std::ostream&                                      os,
    const illumina::Program&                           prog,
    const boost::program_options::options_description& visible,
    const char*                                        msg = nullptr)
{
  usage(
      os,
      prog,
      visible,
      "merge multiple sv-calling statistics files from 'GetAlignmentStats'",
      " > stats",
      msg);
}

void parseMergeAlignmentStatsOptions(
    const illumina::Program& prog, int argc, char* argv[], MergeAlignmentStatsOptions& opt)
{
  namespace po = boost::program_options;
  po::options_description req("configuration");
  // clang-format off
  req.add_options()
  ("align-stats-file", po::value<files_t>(),
   "stats output of 'GetAlignmentStats' (may be specified multiple times)")
  ("output-file", po::value(&opt.outputFilename),
   "write merged stats to filename (default: stdout)");
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
    log_os << "\n" << prog.name() << ": get statistics for SV-calling from alignment files\n\n";
    log_os << "version: " << prog.version() << "\n\n";
    log_os << "usage: " << prog.name() << " [options] > stats\n\n";
    log_os << visible << "\n";
    exit(EXIT_FAILURE);
  }

  std::string errorMsg;
  {
    {
      files_t statsInput;
      if (vm.count("align-stats-file")) {
        statsInput = (boost::any_cast<files_t>(vm["align-stats-file"].value()));
      }
      opt.statsFiles = statsInput;
    }

    errorMsg.clear();
    if (opt.statsFiles.empty()) {
      errorMsg = "Must specify at least one input stats file";
    } else {
      // check that stats files exist, and names do not repeat
      std::set<std::string> nameCheck;
      for (std::string& afile : opt.statsFiles) {
        if (checkAndStandardizeRequiredInputFilePath(afile, "alignment stats file", errorMsg)) break;
        if (nameCheck.count(afile)) {
          std::ostringstream oss;
          oss << "Repeated alignment stats filename: " << afile << "\n";
          errorMsg = oss.str();
          break;
        }
        nameCheck.insert(afile);
      }
    }
  }

  if (!errorMsg.empty()) {
    usage(log_os, prog, visible, errorMsg.c_str());
  }
}
