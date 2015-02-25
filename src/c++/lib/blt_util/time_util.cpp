// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Manta
// Copyright (c) 2013-2014 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/sequencing/licenses/>
//

///
/// \author Chris Saunders
///

#include "io_util.hh"
#include "time_util.hh"

#include <iomanip>
#include <iostream>



void
CpuTimes::
report(
    const double factor,
    const char* tlabel,
    std::ostream& os) const
{
    StreamScoper scoper(os);
    os << std::fixed << std::setprecision(4);
    const double fwall(wall*factor);
    const double fuser(user*factor);
    const double fsystem(system*factor);
    const double total(fuser+fsystem);
    const double perc(100*total/fwall);
    os << fwall << tlabel << " wall, "
       << fuser << tlabel << " user + "
       << fsystem << tlabel << " system = "
       << total << tlabel
       << " CPU (" << std::setprecision(2) << perc << "%)";
}
