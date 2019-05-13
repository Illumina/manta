//
// Manta - Structural Variant and Indel Caller
// Copyright (c) 2013-2019 Illumina, Inc.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//

#include "TestAssemblerOptions.hpp"

#include "blt_util/log.hpp"
#include "common/ProgramUtil.hpp"
#include "options/AlignmentFileOptionsParser.hpp"
#include "options/optionsUtil.hpp"

#include "boost/program_options.hpp"

#include <iostream>

static void usage(
    std::ostream&                                      os,
    const illumina::Program&                           prog,
    const boost::program_options::options_description& visible,
    const char*                                        msg = nullptr)
{
  usage(os, prog, visible, "test manta assembler from command-line", " > contigs", msg);
}

/// \brief Parse TestAssemblerOptions
///
/// \param[out] errorMsg If an error occurs this is set to an end-user targeted error message. Any string
/// content on input is cleared
///
/// \return True if an error occurs while parsing options
static bool parseOptions(
    const boost::program_options::variables_map& vm, TestAssemblerOptions& opt, std::string& errorMsg)
{
  if (parseOptions(vm, opt.alignFileOpt, errorMsg)) return true;
  if (checkAndStandardizeRequiredInputFilePath(opt.referenceFilename, "reference fasta", errorMsg))
    return true;
  return false;
}

void parseTestAssemblerOptions(
    const illumina::Program& prog, int argc, char* argv[], TestAssemblerOptions& opt)
{
  namespace po = boost::program_options;
  po::options_description req("configuration");
  // clang-format off
  req.add_options()
  ("output-file", po::value(&opt.outputFilename),
   "write assembled contigs to filename (default: stdout)")
  ("ref", po::value(&opt.referenceFilename),
   "fasta reference sequence (required)")
  ;
  // clang-format on

  po::options_description help("help");
  help.add_options()("help,h", "print this message");

  po::options_description aligndesc(getOptionsDescription(opt.alignFileOpt));

  po::options_description visible("options");
  visible.add(aligndesc).add(req).add(help);

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
    log_os << "\n" << prog.name() << ": test manta assembler from command-line\n\n";
    log_os << "version: " << prog.version() << "\n\n";
    log_os << "usage: " << prog.name() << " [options] > contigs\n\n";
    log_os << visible << "\n";
    exit(EXIT_FAILURE);
  }

  std::string errorMsg;
  if (parseOptions(vm, opt, errorMsg)) {
    usage(log_os, prog, visible, errorMsg.c_str());
  }
}
