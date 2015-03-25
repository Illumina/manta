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

/// \file

/// \author Chris Saunders
///

#include "blt_util/log.hh"
#include "blt_util/prob_util.hh"

#include <cstdlib>

#include <iomanip>
#include <iostream>



void
check_ln_distro_invalid_value(const char* label,
                              const double val,
                              const unsigned n)
{
    log_os << std::setprecision(14) << std::fixed;
    log_os << "ERROR: " << label << " element [" << n << "] has invalid value: '" << val << "'\n";
    log_os.unsetf(std::ios::fixed);
    exit(EXIT_FAILURE);
}



void
check_ln_distro_invalid_sum(const char* label,
                            const double sum)
{
    log_os << std::setprecision(14) << std::fixed;
    log_os << "ERROR: " << label << " sum is: '" << sum << "'\n";
    log_os.unsetf(std::ios::fixed);
    exit(EXIT_FAILURE);
}
