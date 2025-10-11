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
/// \author Ole Schulz-Trieglaff
///

#pragma once

#include <string>
#include <vector>

/// Information added to each read in the process of assembly
///
struct AssemblyReadInfo {
  bool isUsed = false;

  /// If true, the read was 'used' but filtered out, so there is no meaningful contig id association
  bool isFiltered = false;

  /// If true, the read was an assembled contig
  bool isPseudo = false;

  /// Index of the contigs that this read is used in
  std::vector<unsigned> contigIds;
};

typedef std::vector<std::string>      AssemblyReadInput;
typedef std::vector<bool>             AssemblyReadReversal;
typedef std::vector<AssemblyReadInfo> AssemblyReadOutput;
