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

#include "htsapi/bam_record.hpp"
#include "svgraph/GenomeInterval.hpp"

/// \brief Shared interface for methods which process alignment file regions
///
/// This enables specification of different methods which
/// must traverse a range of reads in a bam file. By abstracting
/// multiple methods to this interface, we can accomplish multiple
/// tasks over a single pass of the BAM records while maintaining
/// isolation of methods
///
struct BamRegionProcessor {
  virtual ~BamRegionProcessor() {}

  /// Provide the index of the next bam file, must be called before switching files/samples
  ///
  /// For each bam index, return the requested interval for this operation,
  /// operations with closely related intervals will be combined
  /// and the union of intervals will be processed.
  virtual const GenomeInterval& nextBamIndex(const unsigned bamIndex) = 0;
};
