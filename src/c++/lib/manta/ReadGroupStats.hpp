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
/// \author Bret Barnes, Xiaoyu Chen
///

#pragma once

#include "blt_util/SizeDistribution.hpp"
#include "common/ReadPairOrient.hpp"
#include "manta/ReadCounter.hpp"

/// Read pair insert stats can be computed for each sample or read group, this
/// class represents the statistics for one group:
///
struct ReadGroupStats {
private:
  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive& ar, const unsigned /*version*/)
  {
    ar& boost::serialization::make_nvp("fragmentSizeDistribution", fragStats);
    ar& boost::serialization::make_nvp("pairOrientation", relOrients);
    ar& boost::serialization::make_nvp("readCount", readCounter);
  }

  ///////////////////////////// data:
public:
  SizeDistribution fragStats;
  ReadPairOrient   relOrients;
  ReadCounter      readCounter;
};

BOOST_CLASS_IMPLEMENTATION(ReadGroupStats, boost::serialization::object_serializable)
