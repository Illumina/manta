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

/// \file

/// \author Chris Saunders
///
#include "blt_util/log.hh"

#include <iostream>

std::ostream& log_os(std::cerr);

void warnOnce(const std::string &msg)
{
    static bool once(false);

    if (!once)
    {
        once = true;
        log_os << msg << std::endl;
    }
}
