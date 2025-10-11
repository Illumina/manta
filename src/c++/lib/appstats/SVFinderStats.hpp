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
