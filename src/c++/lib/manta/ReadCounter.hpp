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
/// \author Naoki Nariai
///

#pragma once

#include <iostream>
#include "boost/foreach.hpp"
#include "boost/serialization/access.hpp"
#include "boost/serialization/nvp.hpp"
#include "boost/serialization/split_member.hpp"

/// \brief Accumulate read statistics scanned for insert size estimation
struct ReadCounter {
  ReadCounter()
    : _totalReadCount(0),
      _totalPairedReadCount(0),
      _totalUnpairedReadCount(0),
      _totalPairedLowMapqReadCount(0),
      _totalHighConfidenceReadPairCount(0)
  {
  }

  unsigned totalReadCount() const { return _totalReadCount; }

  unsigned totalPairedReadCount() const { return _totalPairedReadCount; }

  unsigned totalUnpairedReadCount() const { return _totalUnpairedReadCount; }

  unsigned totalPairedLowMapqReadCount() const { return _totalPairedLowMapqReadCount; }

  unsigned totalHighConfidenceReadPairCount() const { return _totalHighConfidenceReadPairCount; }

  void addReadCount() { _totalReadCount++; }

  void addPairedReadCount() { _totalPairedReadCount++; }

  void addUnpairedReadCount() { _totalUnpairedReadCount++; }

  void addPairedLowMapqReadCount() { _totalPairedLowMapqReadCount++; }

  void addHighConfidenceReadPairCount() { _totalHighConfidenceReadPairCount += 1; }

private:
  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive& ar, const unsigned /*version*/)
  {
    ar& boost::serialization::make_nvp("totalReadCount", _totalReadCount);
    ar& boost::serialization::make_nvp("totalPairedReadCount", _totalPairedReadCount);
    ar& boost::serialization::make_nvp("totalUnpairedReadCount", _totalUnpairedReadCount);
    ar& boost::serialization::make_nvp("totalPairedLowMapqReadCount", _totalPairedLowMapqReadCount);
    ar& boost::serialization::make_nvp("totalHighConfidenceReadPairCount", _totalHighConfidenceReadPairCount);
  }

  ///////////////////////////////////// data:
  unsigned _totalReadCount;
  unsigned _totalPairedReadCount;
  unsigned _totalUnpairedReadCount;
  unsigned _totalPairedLowMapqReadCount;
  unsigned _totalHighConfidenceReadPairCount;
};

BOOST_CLASS_IMPLEMENTATION(ReadCounter, object_serializable)

inline std::ostream& operator<<(std::ostream& os, const ReadCounter& rs)
{
  os << "\tTotal sampled reads: " + std::to_string(rs.totalReadCount()) + "\n"
     << "\tTotal sampled paired reads: " + std::to_string(rs.totalPairedReadCount()) + "\n"
     << "\tTotal sampled paired reads passing MAPQ filter: " +
            std::to_string(rs.totalPairedReadCount() - rs.totalPairedLowMapqReadCount()) + "\n"
     << "\tTotal sampled high-confidence read pairs passing all filters: " +
            std::to_string(rs.totalHighConfidenceReadPairCount()) + "\n";
  return os;
}
