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

/// Input parameters for SmallAssembler, a simple local de-bruijn graph assembler
///
struct SmallAssemblerOptions {
  SmallAssemblerOptions() {}

  /// the symbol set used during assembly
  std::string alphabet = "ACGT";

  /// minimum basecall quality for assembly input
  uint8_t minQval = 5;

  /// initial word (kmer) length
  unsigned minWordLength   = 41;
  unsigned maxWordLength   = 76;
  unsigned wordStepSize    = 5;
  unsigned minContigLength = 15;

  /// min. coverage required for contig extension
  unsigned minCoverage = 1;

  /// coverage required for conservative contig sub-range
  unsigned minConservativeCoverage = 2;

  /// max error rates allowed during contig extension
  double maxError = 0.35;

  /// min. number of reads required to start assembly
  unsigned minSeedReads = 3;

  /// Max. number of assembly iterations for a cluster before we give up
  unsigned maxAssemblyIterations = 10;
};
