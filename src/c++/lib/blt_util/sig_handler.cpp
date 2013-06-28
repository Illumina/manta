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

/// \author Chris Saunders
///

#include "blt_util/log.hh"
#include "blt_util/sig_handler.hh"

#include <cstdlib>
#include <signal.h>

#include <iostream>
#include <string>


static std::string _progname;
static std::string _cmdline;


static
void
blt_sig_handler (int sig)
{

    switch (sig)
    {
    case SIGTERM:
        log_os << "ERROR: " << _progname << " received termination signal. cmdline: " << _cmdline << std::endl;
        exit(EXIT_FAILURE);
#ifndef _WIN32
    case SIGINT:
        log_os << "ERROR: " << _progname << " received interupt signal. cmdline: " << _cmdline << std::endl;
        exit(EXIT_FAILURE);
#endif
    default:
        log_os << "INFO: " << _progname << " received signal no: " << sig << std::endl;
        break;
    }
}



void
initialize_blt_signals(const char* progname,
                       const char* cmdline)
{

    _progname=progname;
    _cmdline=cmdline;

    signal(SIGTERM, blt_sig_handler);
#ifndef _WIN32
    signal(SIGINT, blt_sig_handler);
#endif
}
