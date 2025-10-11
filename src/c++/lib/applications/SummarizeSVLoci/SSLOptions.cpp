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

#include "SSLOptions.hpp"

#include "blt_util/log.hpp"
#include "boost/filesystem.hpp"
#include "boost/program_options.hpp"

#include <iostream>
#include "../../common/ProgramUtil.hpp"

static void usage(
    std::ostream&                                      os,
    const illumina::Program&                           prog,
    const boost::program_options::options_description& visible,
    const char*                                        msg = nullptr)
{
  usage(os, prog, visible, "write graph summary stats to stdout (tsv format)", " [ > output ]", msg);
}

void parseSSLOptions(const illumina::Program& prog, int argc, char* argv[], SSLOptions& opt)
{
  namespace po = boost::program_options;
  po::options_description req("configuration");
  // clang-format off
  req.add_options()
  ("graph-file", po::value(&opt.graphFilename),
   "sv locus graph file")
  ("global", po::value(&opt.isGlobalStats)->zero_tokens(),
   "provide global stats on full graph (default output is per-locus stats)")
  ("output-file", po::value(&opt.outputFilename),
   "write graph summary stats to filename (default: stdout)");
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

  // fast check of config state:
  if (opt.graphFilename.empty()) {
    usage(log_os, prog, visible, "Must specify sv locus graph file");
  }
  if (!boost::filesystem::exists(opt.graphFilename)) {
    usage(log_os, prog, visible, "SV locus graph file does not exist");
  }
}
