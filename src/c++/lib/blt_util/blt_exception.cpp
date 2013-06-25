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

/// \file

/// \author Chris Saunders
///
#include "blt_util/blt_exception.hh"

#ifdef KILL_EXCEPTIONS
#include "blt_util/log.hh"

#include <cstdlib>

#include <iostream>
#endif



blt_exception::
blt_exception(const char* s)
    : message(s)
{
#ifdef KILL_EXCEPTIONS
    log_os << "ERROR:: " << s << std::endl;
    abort();
#endif
}

