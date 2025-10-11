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

#include "blt_util/prob_util.hpp"
#include "blt_util/log.hpp"

#include <cstdlib>

#include <iomanip>
#include <iostream>

void check_ln_distro_invalid_value(const char* label, const double val, const unsigned n)
{
  log_os << std::setprecision(14) << std::fixed;
  log_os << "ERROR: " << label << " element [" << n << "] has invalid value: '" << val << "'\n";
  log_os.unsetf(std::ios::fixed);
  exit(EXIT_FAILURE);
}

void check_ln_distro_invalid_sum(const char* label, const double sum)
{
  log_os << std::setprecision(14) << std::fixed;
  log_os << "ERROR: " << label << " sum is: '" << sum << "'\n";
  log_os.unsetf(std::ios::fixed);
  exit(EXIT_FAILURE);
}
