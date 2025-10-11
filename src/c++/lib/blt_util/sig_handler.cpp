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

#include "blt_util/sig_handler.hpp"
#include "blt_util/log.hpp"

#include <signal.h>
#include <cstdlib>

#include <iostream>
#include <string>

static std::string _progname;
static std::string _cmdline;

static void blt_sig_handler(int sig)
{
  switch (sig) {
  case SIGTERM:
    log_os << "ERROR: " << _progname << " received termination signal. cmdline: " << _cmdline << std::endl;
    exit(EXIT_FAILURE);
#ifndef _WIN32
  case SIGINT:
    log_os << "ERROR: " << _progname << " received interrupt signal. cmdline: " << _cmdline << std::endl;
    exit(EXIT_FAILURE);
#endif
  default:
    log_os << "INFO: " << _progname << " received signal no: " << sig << std::endl;
    break;
  }
}

void initialize_blt_signals(const char* progname, const char* cmdline)
{
  _progname = progname;
  _cmdline  = cmdline;

  signal(SIGTERM, blt_sig_handler);
#ifndef _WIN32
  signal(SIGINT, blt_sig_handler);
#endif
}
