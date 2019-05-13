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
