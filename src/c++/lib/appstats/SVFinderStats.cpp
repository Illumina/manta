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

#include "SVFinderStats.hpp"
#include <iostream>

void SVFinderStats::report(std::ostream& os) const
{
  os << "EdgeFilter\t" << edgeFilter << "\n";
  os << "SemiMappedFilter\t" << semiMappedFilter << "\n";
  os << "ComplexLowCountFilter\t" << ComplexLowCountFilter << "\n";
  os << "ComplexLowSignalFilter\t" << ComplexLowSignalFilter << "\n";
  os << "UnmatchedReadPairFilter\t" << unmatchedReadPairFilter << "\n";
}
