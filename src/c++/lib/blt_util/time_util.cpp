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

#include "time_util.hpp"
#include "io_util.hpp"

#include <iomanip>
#include <iostream>

void CpuTimes::report(const double factor, const char* tlabel, std::ostream& os) const
{
  StreamScoper scoper(os);
  os << std::fixed << std::setprecision(4);
  const double fwall(wall * factor);
  const double fuser(user * factor);
  const double fsystem(system * factor);
  const double total(fuser + fsystem);
  const double perc(100 * total / fwall);
  os << fwall << tlabel << " wall, " << fuser << tlabel << " user + " << fsystem << tlabel
     << " system = " << total << tlabel << " CPU (" << std::setprecision(2) << perc << "%)";
}
