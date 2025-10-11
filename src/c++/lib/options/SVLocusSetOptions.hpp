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

#pragma once

/// parameters specific to SVLocusSet:
///
struct SVLocusSetOptions {
  explicit SVLocusSetOptions(const unsigned initObservationWeight = 1)
    : observationWeight(initObservationWeight),
      minMergeEdgeObservations(3),
      maxSearchCount(500),
      maxSearchDensity(0.5)
  {
  }

  unsigned getMinMergeEdgeCount() const { return (observationWeight * minMergeEdgeObservations); }

  template <class Archive>
  void serialize(Archive& ar, const unsigned /* version */)
  {
    ar& observationWeight;
    ar& minMergeEdgeObservations;
    ar& maxSearchCount;
    ar& maxSearchDensity;
  }

  /// Used to translate graph edges counts to observations
  unsigned observationWeight;

  /// To reduce noise in the graph, we only merge once shared edges reach this number of observations
  unsigned minMergeEdgeObservations;

  /// The search for intersecting regions in the graph stops once this number is reached
  unsigned maxSearchCount;

  /// The search for intersecting regions in the graph stops once this many regions/base are found
  float maxSearchDensity;
};
