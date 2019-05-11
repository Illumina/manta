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

/// \file
/// \author Chris Saunders
///

#pragma once

#include "boost/serialization/nvp.hpp"

#include <cstdint>

#include <iosfwd>

/// statistics accumulated during the SVFinder process
///
/// this includes trackers for filtered edges, sv candidates, and reads/read-pairs
///
struct SVFinderStats {
  SVFinderStats() {}

  void merge(const SVFinderStats& rhs)
  {
    edgeFilter += rhs.edgeFilter;
    semiMappedFilter += rhs.semiMappedFilter;
    ComplexLowCountFilter += rhs.ComplexLowCountFilter;
    ComplexLowSignalFilter += rhs.ComplexLowSignalFilter;
    unmatchedReadPairFilter += rhs.unmatchedReadPairFilter;
  }

  template <class Archive>
  void serialize(Archive& ar, const unsigned /* version */)
  {
    ar& BOOST_SERIALIZATION_NVP(edgeFilter) & BOOST_SERIALIZATION_NVP(semiMappedFilter) &
        BOOST_SERIALIZATION_NVP(ComplexLowCountFilter) & BOOST_SERIALIZATION_NVP(ComplexLowSignalFilter) &
        BOOST_SERIALIZATION_NVP(unmatchedReadPairFilter);
  }

  void report(std::ostream& os) const;

  /// number of edges filtered out from the SV finding process
  uint64_t edgeFilter = 0;

  /// number of sv candidates filtered out for being supported only be semi-mapped read pairs
  uint64_t semiMappedFilter = 0;

  /// number of complex candidates filtered out for not meeting hard evidence count thresholds
  uint64_t ComplexLowCountFilter = 0;

  /// number of complex candidates filtered out for evidence counts which are not significant relative to
  /// noise rates in the data
  uint64_t ComplexLowSignalFilter = 0;

  /// number of read-pairs filtered out becuase the read1 and read2 alignment records contained conflicting
  /// information
  uint64_t unmatchedReadPairFilter = 0;
};

BOOST_CLASS_IMPLEMENTATION(SVFinderStats, boost::serialization::object_serializable)
